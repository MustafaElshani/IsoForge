// In-process spliced alignment with minimap2.
//
// Reads (FASTA or FASTQ, gzipped or not) are mapped with a minimap2 preset, and each record is the SAM text minimap2's
// writer (mm_write_sam3) produces, parsed by htslib. The index is built once (or loaded from a saved .mmi) and serves
// every read file.

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int};
use std::path::{Path, PathBuf};
use std::ptr;

use minimap2_sys as mm;
use rayon::prelude::*;
use rust_htslib::bam;

use crate::errors::{IsoForgeError, IsoForgeResult};

/// How reads are aligned; each field is the minimap2 option of the same meaning.
#[derive(Debug, Clone)]
pub struct AlignOptions {
    /// `-x`: `splice` (noisy long reads) or `splice:hq` (accurate long reads).
    pub preset: String,
    /// `-k`: minimizer length (the preset's when None).
    pub k: Option<i32>,
    /// `-uf`: match splice signals on the forward strand only, for reads already in transcript orientation;
    /// `-ub` (both strands) otherwise.
    pub forward_only: bool,
    /// Print secondary alignments (`--secondary=yes`).
    pub secondary: bool,
    /// `--MD`.
    pub md: bool,
    /// `--spsc`: splice-site scores, and `--spsc-scale`.
    pub splice_scores: Option<PathBuf>,
    pub splice_score_scale: f32,
    /// `--junc-bed`: annotated junctions (BED12 or minimap2's junction BED).
    pub junc_bed: Option<PathBuf>,
    /// `-y`: copy the FASTA/FASTQ comment to the SAM record.
    pub copy_comment: bool,
    /// `-d`: where to save the index built from a FASTA (a path ending .mmi is read as a saved index instead).
    pub index_out: Option<PathBuf>,
}

impl Default for AlignOptions {
    fn default() -> Self {
        Self { preset: "splice".into(), k: None, forward_only: false, secondary: false, md: true, splice_scores: None,
               splice_score_scale: 0.7, junc_bed: None, copy_comment: false, index_out: None }
    }
}

/// Counts from aligning one file.
#[derive(Debug, Clone, Copy, Default)]
pub struct AlignStats {
    pub reads: u64,
    pub unmapped: u64,
    pub primary: u64,
    pub supplementary: u64,
    pub secondary: u64,
}

/// A minimap2 index with its mapping options, ready to align reads from any number of files.
pub struct Aligner {
    mi: *mut mm::mm_idx_t,
    opt: mm::mm_mapopt_t,
    header: bam::HeaderView,
    header_text: String,
}

// The index and options are only read while mapping; every thread maps with its own buffer.
unsafe impl Send for Aligner {}
unsafe impl Sync for Aligner {}

impl Drop for Aligner {
    fn drop(&mut self) {
        if !self.mi.is_null() { unsafe { mm::mm_idx_destroy(self.mi) } }
    }
}

fn cstring(p: &Path) -> IsoForgeResult<CString> {
    CString::new(p.to_string_lossy().as_bytes()).map_err(|e| IsoForgeError::ConfigError(e.to_string()))
}

/// minimap2's read records, shared read-only across the mapping threads.
struct Batch(*const mm::mm_bseq1_t, usize);
unsafe impl Send for Batch {}
unsafe impl Sync for Batch {}

impl Batch {
    fn get(&self, i: usize) -> &mm::mm_bseq1_t { assert!(i < self.1); unsafe { &*self.0.add(i) } }
}

/// One thread's minimap2 buffer (freed when dropped).
struct ThreadBuffer(*mut mm::mm_tbuf_t, mm::kstring_t);

impl ThreadBuffer {
    fn new() -> Self { Self(unsafe { mm::mm_tbuf_init() }, mm::kstring_t { l: 0, m: 0, s: ptr::null_mut() }) }
}

impl Drop for ThreadBuffer {
    fn drop(&mut self) {
        unsafe {
            mm::mm_tbuf_destroy(self.0);
            if !self.1.s.is_null() { libc_free(self.1.s as *mut std::ffi::c_void); }
        }
    }
}

extern "C" {
    #[link_name = "free"]
    fn libc_free(p: *mut std::ffi::c_void);
    /// Loads splice-site scores into a minimap2 index, scaled by `scale`.
    fn mm_idx_spsc_read2(idx: *mut mm::mm_idx_t, fn_: *const c_char, max_sc: i32, scale: f32) -> i32;
}

/// A minimap2 index of a few sequences held in memory, for asking where another sequence lies in them: `cover` maps a
/// query with base-level alignment and reports the query bases its primary hits cover and the matching bases in them.
pub struct SequenceIndex {
    mi: *mut mm::mm_idx_t,
    opt: mm::mm_mapopt_t,
}

unsafe impl Send for SequenceIndex {}
unsafe impl Sync for SequenceIndex {}

impl Drop for SequenceIndex {
    fn drop(&mut self) {
        if !self.mi.is_null() { unsafe { mm::mm_idx_destroy(self.mi) } }
    }
}

impl SequenceIndex {
    /// The index of `targets` with minimap2's `preset`. None when no target holds sequence.
    pub fn new(targets: &[Vec<u8>], preset: &str) -> Option<Self> {
        let seqs: Vec<CString> = targets.iter().filter(|t| !t.is_empty())
            .filter_map(|t| CString::new(t.iter().map(|&b| if b == 0 { b'N' } else { b }).collect::<Vec<u8>>()).ok()).collect();
        if seqs.is_empty() { return None; }
        let names: Vec<CString> = (0..seqs.len()).map(|i| CString::new(format!("t{i}")).unwrap()).collect();
        unsafe {
            let mut ipt: mm::mm_idxopt_t = std::mem::zeroed();
            let mut opt: mm::mm_mapopt_t = std::mem::zeroed();
            mm::mm_set_opt(ptr::null(), &mut ipt, &mut opt);
            let preset = CString::new(preset).ok()?;
            if mm::mm_set_opt(preset.as_ptr(), &mut ipt, &mut opt) < 0 { return None; }
            opt.flag |= mm::MM_F_CIGAR as i64;
            let mut sp: Vec<*const c_char> = seqs.iter().map(|c| c.as_ptr()).collect();
            let mut np: Vec<*const c_char> = names.iter().map(|c| c.as_ptr()).collect();
            let mi = mm::mm_idx_str(ipt.w as c_int, ipt.k as c_int, (ipt.flag & 1) as c_int, ipt.bucket_bits as c_int,
                                    sp.len() as c_int, sp.as_mut_ptr(), np.as_mut_ptr());
            if mi.is_null() { return None; }
            mm::mm_mapopt_update(&mut opt, mi);
            Some(Self { mi, opt })
        }
    }

    /// (query bases covered by the primary hits, matching bases in them) for `query`.
    pub fn cover(&self, query: &[u8]) -> (usize, usize) {
        if query.is_empty() { return (0, 0); }
        let Ok(q) = CString::new(query.iter().map(|&b| if b == 0 { b'N' } else { b }).collect::<Vec<u8>>()) else { return (0, 0) };
        let name = CString::new("q").unwrap();
        let mut spans: Vec<(i32, i32)> = Vec::new();
        let mut matches = 0usize;
        unsafe {
            let tbuf = mm::mm_tbuf_init();
            let mut n: c_int = 0;
            let regs = mm::mm_map(self.mi, query.len() as c_int, q.as_ptr(), &mut n, tbuf, &self.opt, name.as_ptr());
            for i in 0..n.max(0) as usize {
                let r = &*regs.add(i);
                if r.id == r.parent {
                    spans.push((r.qs, r.qe));
                    matches += r.mlen.max(0) as usize;
                }
                if !r.p.is_null() { libc_free(r.p as *mut std::ffi::c_void); }
            }
            if !regs.is_null() { libc_free(regs as *mut std::ffi::c_void); }
            mm::mm_tbuf_destroy(tbuf);
        }
        spans.sort_unstable();
        let (mut covered, mut end) = (0usize, i32::MIN);
        for (s, e) in spans {
            let s = s.max(end);
            if e > s { covered += (e - s) as usize; end = e; }
        }
        (covered, matches)
    }
}

impl Aligner {
    /// Builds the index of `genome` (FASTA, or a saved .mmi) with `opts` and loads the annotated junctions and
    /// splice-site scores into it.
    pub fn new(genome: &Path, opts: &AlignOptions, threads: usize) -> IsoForgeResult<Self> {
        unsafe {
            mm::mm_verbose = 1;
            let mut ipt: mm::mm_idxopt_t = std::mem::zeroed();
            let mut opt: mm::mm_mapopt_t = std::mem::zeroed();
            mm::mm_set_opt(ptr::null(), &mut ipt, &mut opt);
            let preset = CString::new(opts.preset.as_str()).map_err(|e| IsoForgeError::ConfigError(e.to_string()))?;
            if mm::mm_set_opt(preset.as_ptr(), &mut ipt, &mut opt) < 0 {
                return Err(IsoForgeError::ConfigError(format!("unknown minimap2 preset {:?}", opts.preset)));
            }
            opt.flag |= (mm::MM_F_OUT_SAM | mm::MM_F_CIGAR) as i64;
            if let Some(k) = opts.k { ipt.k = k as _; }
            if opts.md { opt.flag |= mm::MM_F_OUT_MD as i64; }
            if opts.forward_only {
                opt.flag |= mm::MM_F_SPLICE_FOR as i64;
                opt.flag &= !(mm::MM_F_SPLICE_REV as i64);
            }
            if !opts.secondary { opt.flag |= mm::MM_F_NO_PRINT_2ND as i64; }
            if opts.copy_comment { opt.flag |= mm::MM_F_COPY_COMMENT as i64; }
            if mm::mm_check_opt(&ipt, &opt) < 0 {
                return Err(IsoForgeError::ConfigError("inconsistent minimap2 options".into()));
            }
            let fa = cstring(genome)?;
            let out = match &opts.index_out { Some(p) => Some(cstring(p)?), None => None };
            let rdr = mm::mm_idx_reader_open(fa.as_ptr(), &ipt, out.as_ref().map_or(ptr::null(), |c| c.as_ptr()));
            if rdr.is_null() {
                return Err(IsoForgeError::ConfigError(format!("cannot read the genome {}", genome.display())));
            }
            let mi = mm::mm_idx_reader_read(rdr, threads.max(1) as c_int);
            let single = mm::mm_idx_reader_eof(rdr) != 0;
            mm::mm_idx_reader_close(rdr);
            if mi.is_null() {
                return Err(IsoForgeError::ConfigError(format!("no sequences in {}", genome.display())));
            }
            if !single {
                mm::mm_idx_destroy(mi);
                return Err(IsoForgeError::ConfigError("the genome index is split into parts; raise the index batch size".into()));
            }
            mm::mm_mapopt_update(&mut opt, mi);
            if let Some(bed) = &opts.junc_bed {
                let c = cstring(bed)?;
                mm::mm_idx_bed_read(mi, c.as_ptr(), 1);
                if (*mi).I.is_null() {
                    mm::mm_idx_destroy(mi);
                    return Err(IsoForgeError::ConfigError(format!("cannot read the junction BED {}", bed.display())));
                }
            }
            if let Some(sc) = &opts.splice_scores {
                let c = cstring(sc)?;
                mm_idx_spsc_read2(mi, c.as_ptr(), mm::mm_max_spsc_bonus(&opt), opts.splice_score_scale);
                if (*mi).spsc.is_null() {
                    mm::mm_idx_destroy(mi);
                    return Err(IsoForgeError::ConfigError(format!("cannot read the splice scores {}", sc.display())));
                }
            }
            let mut header_text = String::from("@HD\tVN:1.6\tSO:unsorted\n");
            for i in 0..(*mi).n_seq as usize {
                let s = &*(*mi).seq.add(i);
                header_text.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", CStr::from_ptr(s.name).to_string_lossy(), s.len));
            }
            header_text.push_str(&format!("@PG\tID:minimap2\tPN:minimap2\tVN:2.30-r1287\tCL:isoforge align -x {}{}{}{}{}\n",
                opts.preset, opts.k.map_or(String::new(), |k| format!(" -k{k}")), if opts.forward_only { " -uf" } else { "" },
                if opts.md { " --MD" } else { "" }, if opts.splice_scores.is_some() { " --spsc" } else { "" }));
            let header = bam::HeaderView::from_bytes(header_text.as_bytes());
            // builds the header's name lookup, before records are parsed in parallel
            if (*mi).n_seq > 0 {
                let first = CStr::from_ptr((*(*mi).seq).name).to_bytes().to_vec();
                header.tid(&first);
            }
            Ok(Self { mi, opt, header, header_text })
        }
    }

    /// The BAM header of the records (the genome's sequences, in index order).
    pub fn header(&self) -> &bam::HeaderView { &self.header }
    pub fn header_text(&self) -> &str { &self.header_text }

    /// The SAM lines of one read, as minimap2 writes them.
    unsafe fn map_one(&self, t: &mm::mm_bseq1_t, buf: &mut ThreadBuffer, out: &mut Vec<u8>) -> (u64, u64, u64) {
        let mut n_reg: c_int = 0;
        let mut regs: *mut mm::mm_reg1_t = ptr::null_mut();
        let qlen: c_int = t.l_seq;
        let qseq: *const c_char = t.seq;
        mm::mm_map_frag(self.mi, 1, &qlen, &qseq as *const *const c_char as *mut *const c_char, &mut n_reg, &mut regs,
                        buf.0, &self.opt, t.name);
        let rep_len = (*buf.0).rep_len;
        let regss: *const *const mm::mm_reg1_t = &(regs as *const mm::mm_reg1_t);
        let (mut prim, mut supp, mut sec) = (0u64, 0u64, 0u64);
        let mut write = |j: c_int, out: &mut Vec<u8>| {
            mm::mm_write_sam3(&mut buf.1, self.mi, t, 0, j, 1, &n_reg, regss, ptr::null_mut(), self.opt.flag, rep_len);
            out.extend_from_slice(std::slice::from_raw_parts(buf.1.s as *const u8, buf.1.l as usize));
            out.push(b'\n');
        };
        for j in 0..n_reg {
            let r = &*regs.add(j as usize);
            if r.id != r.parent {
                if self.opt.flag & mm::MM_F_NO_PRINT_2ND as i64 != 0 { continue; }
                sec += 1;
            } else if r.sam_pri() != 0 { prim += 1; } else { supp += 1; }
            write(j, out);
        }
        if n_reg == 0 { write(-1, out); }
        for j in 0..n_reg as usize { libc_free((*regs.add(j)).p as *mut std::ffi::c_void); }
        if !regs.is_null() { libc_free(regs as *mut std::ffi::c_void); }
        (prim, supp, sec)
    }

    /// Aligns every read of `reads` and hands the records to `sink` batch by batch, in input order.
    pub fn align_file<F>(&self, reads: &Path, mut sink: F) -> IsoForgeResult<AlignStats>
    where F: FnMut(Vec<bam::Record>) -> IsoForgeResult<()> {
        let path = cstring(reads)?;
        let fp = unsafe { mm::mm_bseq_open(path.as_ptr()) };
        if fp.is_null() {
            return Err(IsoForgeError::ConfigError(format!("cannot open reads {}", reads.display())));
        }
        let with_comment = (self.opt.flag & mm::MM_F_COPY_COMMENT as i64 != 0) as c_int;
        let mut stats = AlignStats::default();
        let result = (|| -> IsoForgeResult<()> {
            loop {
                let mut n: c_int = 0;
                let seqs = unsafe { mm::mm_bseq_read3(fp, self.opt.mini_batch_size, 1, with_comment, 0, &mut n) };
                if seqs.is_null() || n == 0 { break; }
                let shared = Batch(seqs, n as usize);
                // SAM text per read, then records, in input order
                let texts: Vec<(Vec<u8>, (u64, u64, u64))> = (0..shared.1).into_par_iter().with_min_len(16)
                    .map_init(ThreadBuffer::new, |buf, i| {
                        let t = shared.get(i);
                        let mut out = Vec::with_capacity(4 * t.l_seq as usize + 256);
                        let c = unsafe { self.map_one(t, buf, &mut out) };
                        (out, c)
                    }).collect();
                let batch: &[mm::mm_bseq1_t] = unsafe { std::slice::from_raw_parts(seqs, n as usize) };
                let records: Vec<IsoForgeResult<Vec<bam::Record>>> = texts.par_iter().map(|(text, _)| {
                    text.split(|&b| b == b'\n').filter(|l| !l.is_empty()).map(|line| {
                        bam::Record::from_sam(&self.header, line)
                            .map_err(|e| IsoForgeError::HtsError(format!("SAM record from minimap2: {e}")))
                    }).collect()
                }).collect();
                for (_, (p, s, x)) in &texts {
                    stats.reads += 1;
                    if p + s + x == 0 { stats.unmapped += 1; }
                    stats.primary += p; stats.supplementary += s; stats.secondary += x;
                }
                unsafe {
                    for t in batch {
                        for p in [t.name, t.seq, t.qual, t.comment] {
                            if !p.is_null() { libc_free(p as *mut std::ffi::c_void); }
                        }
                    }
                    libc_free(seqs as *mut std::ffi::c_void);
                }
                let mut all = Vec::new();
                for r in records { all.extend(r?); }
                sink(all)?;
            }
            Ok(())
        })();
        unsafe { mm::mm_bseq_close(fp) };
        result.map(|_| stats)
    }

    /// Aligns `reads` into an unsorted BAM at `out`.
    pub fn align_to_bam(&self, reads: &Path, out: &Path, threads: usize) -> IsoForgeResult<AlignStats> {
        let header = bam::Header::from_template(&self.header);
        let mut w = bam::Writer::from_path(out, &header, bam::Format::Bam)
            .map_err(|e| IsoForgeError::HtsError(format!("{}: {e}", out.display())))?;
        w.set_threads(threads.clamp(1, 8)).map_err(|e| IsoForgeError::HtsError(e.to_string()))?;
        self.align_file(reads, |records| {
            for r in &records { w.write(r).map_err(|e| IsoForgeError::HtsError(e.to_string()))?; }
            Ok(())
        })
    }
}

/// The coordinate sort key of a record: reference (unmapped last), position, reverse strand; records tied on all three
/// keep their input order.
fn coord_key(r: &bam::Record) -> (u64, u64, bool) {
    let tid = if r.tid() < 0 { u64::MAX } else { r.tid() as u64 };
    let pos = if r.pos() < 0 { u64::MAX } else { r.pos() as u64 };
    (tid, pos, r.is_reverse())
}

/// Bytes a record holds in memory, for the sort's budget.
fn record_bytes(r: &bam::Record) -> usize { r.inner().l_data.max(0) as usize + 96 }

impl Aligner {
    /// Aligns `reads` into a coordinate-sorted BAM at `out`, with its index. Records are sorted `mem_bytes` at a time;
    /// each full batch is written to `tmp_dir` (fast compression) and the batches are merged, ties in batch order.
    pub fn align_to_sorted_bam(&self, reads: &Path, out: &Path, tmp_dir: &Path, threads: usize, mem_bytes: usize)
                               -> IsoForgeResult<AlignStats> {
        let text = self.header_text.replacen("SO:unsorted", "SO:coordinate", 1);
        let header = bam::Header::from_template(&bam::HeaderView::from_bytes(text.as_bytes()));
        let hts = |e: rust_htslib::errors::Error| IsoForgeError::HtsError(e.to_string());
        let write_sorted = |mut recs: Vec<bam::Record>, path: &Path, fast: bool| -> IsoForgeResult<()> {
            recs.par_sort_by_key(coord_key);
            let mut w = bam::Writer::from_path(path, &header, bam::Format::Bam).map_err(hts)?;
            w.set_threads(threads.clamp(1, 16)).map_err(hts)?;
            if fast { w.set_compression_level(bam::CompressionLevel::Fastest).map_err(hts)?; }
            for r in &recs { w.write(r).map_err(hts)?; }
            Ok(())
        };
        std::fs::create_dir_all(tmp_dir)?;
        let stem = out.file_stem().map_or("isoforge".into(), |s| s.to_string_lossy().to_string());
        let mut chunks: Vec<PathBuf> = Vec::new();
        let (mut buf, mut used): (Vec<bam::Record>, usize) = (Vec::new(), 0);
        let stats = self.align_file(reads, |records| {
            for r in records {
                used += record_bytes(&r);
                buf.push(r);
            }
            if used >= mem_bytes {
                let path = tmp_dir.join(format!("{stem}.sort.{:04}.bam", chunks.len()));
                write_sorted(std::mem::take(&mut buf), &path, true)?;
                chunks.push(path);
                used = 0;
            }
            Ok(())
        })?;
        if chunks.is_empty() {
            write_sorted(buf, out, false)?;
        } else {
            if !buf.is_empty() {
                let path = tmp_dir.join(format!("{stem}.sort.{:04}.bam", chunks.len()));
                write_sorted(std::mem::take(&mut buf), &path, true)?;
                chunks.push(path);
            }
            // k-way merge: the smallest key first, ties to the earlier batch
            use bam::Read;
            let mut readers: Vec<bam::Reader> = chunks.iter().map(|p| bam::Reader::from_path(p).map_err(hts)).collect::<Result<_, _>>()?;
            for r in readers.iter_mut() { r.set_threads(2).map_err(hts)?; }
            let mut heads: Vec<Option<bam::Record>> = Vec::with_capacity(readers.len());
            let mut heap = std::collections::BinaryHeap::new();
            for (i, r) in readers.iter_mut().enumerate() {
                let mut rec = bam::Record::new();
                match r.read(&mut rec) {
                    Some(Ok(())) => { heap.push(std::cmp::Reverse((coord_key(&rec), i))); heads.push(Some(rec)); }
                    Some(Err(e)) => return Err(hts(e)),
                    None => heads.push(None),
                }
            }
            let mut w = bam::Writer::from_path(out, &header, bam::Format::Bam).map_err(hts)?;
            w.set_threads(threads.clamp(1, 16)).map_err(hts)?;
            while let Some(std::cmp::Reverse((_, i))) = heap.pop() {
                let rec = heads[i].take().expect("a head for every heap entry");
                w.write(&rec).map_err(hts)?;
                let mut next = rec;
                match readers[i].read(&mut next) {
                    Some(Ok(())) => { heap.push(std::cmp::Reverse((coord_key(&next), i))); heads[i] = Some(next); }
                    Some(Err(e)) => return Err(hts(e)),
                    None => {}
                }
            }
            drop(w);
            for p in &chunks { std::fs::remove_file(p).ok(); }
        }
        bam::index::build(out, None, bam::index::Type::Bai, threads.clamp(1, 16) as u32).map_err(hts)?;
        Ok(stats)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn a_sequence_index_finds_an_exon_copy_and_not_unrelated_sequence() {
        let mut x: u64 = 12345;
        let mut base = || { x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); b"ACGT"[(x >> 33) as usize % 4] };
        let exon: Vec<u8> = (0..180).map(|_| base()).collect();
        let intron: Vec<u8> = (0..3000).map(|_| base()).collect();
        let other: Vec<u8> = (0..180).map(|_| base()).collect();
        // a gene holding a copy of the exon after 3 kb of other sequence
        let mut gene = intron.clone();
        gene.extend_from_slice(&exon);
        gene.extend((0..2000).map(|_| base()));
        let ix = super::SequenceIndex::new(&[gene], "splice:hq").unwrap();
        let (cov, m) = ix.cover(&exon);
        assert!(cov >= 170 && m >= 170, "{cov} {m}");
        assert_eq!(ix.cover(&other).0, 0);
    }

    use super::*;
    use std::io::Write;

    /// A two-exon gene on a random 20 kb contig; a read spliced across its intron maps with an N in its CIGAR.
    #[test]
    fn a_spliced_read_maps_with_its_intron() {
        let dir = std::env::temp_dir().join(format!("isoforge_align_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut state = 12345u64;
        let mut base = || { state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); b"ACGT"[(state >> 33) as usize % 4] };
        let mut contig: Vec<u8> = (0..20_000).map(|_| base()).collect();
        // exon 1 [5000, 5300), intron with GT..AG, exon 2 [6300, 6600)
        contig[5300] = b'G'; contig[5301] = b'T'; contig[6298] = b'A'; contig[6299] = b'G';
        let fa = dir.join("g.fa");
        let mut f = std::fs::File::create(&fa).unwrap();
        writeln!(f, ">chrT").unwrap();
        for chunk in contig.chunks(60) { f.write_all(chunk).unwrap(); f.write_all(b"\n").unwrap(); }
        let read: Vec<u8> = contig[5000..5300].iter().chain(&contig[6300..6600]).copied().collect();
        let fq = dir.join("r.fa");
        std::fs::write(&fq, format!(">r1\n{}\n", String::from_utf8_lossy(&read))).unwrap();
        let aligner = Aligner::new(&fa, &AlignOptions::default(), 1).unwrap();
        let mut got = Vec::new();
        let stats = aligner.align_file(&fq, |r| { got.extend(r); Ok(()) }).unwrap();
        assert_eq!((stats.reads, stats.primary, stats.unmapped), (1, 1, 0));
        let r = &got[0];
        assert_eq!(r.pos(), 5000);
        assert_eq!(r.cigar().to_string(), "300M1000N300M");
        assert!(r.aux(b"MD").is_ok());
        // sorted: two reads out of order, spilled one batch at a time (a 1-byte budget), come out by position
        let read2: Vec<u8> = contig[12000..12400].to_vec();
        std::fs::write(&fq, format!(">r2\n{}\n>r1\n{}\n", String::from_utf8_lossy(&read2), String::from_utf8_lossy(&read))).unwrap();
        let out = dir.join("s.bam");
        aligner.align_to_sorted_bam(&fq, &out, &dir.join("tmp"), 1, 1).unwrap();
        use rust_htslib::bam::Read;
        let names: Vec<String> = bam::Reader::from_path(&out).unwrap().records()
            .map(|r| String::from_utf8_lossy(r.unwrap().qname()).to_string()).collect();
        assert_eq!(names, vec!["r1", "r2"]);
        assert!(dir.join("s.bam.bai").exists());
        std::fs::remove_dir_all(&dir).ok();
    }
}
