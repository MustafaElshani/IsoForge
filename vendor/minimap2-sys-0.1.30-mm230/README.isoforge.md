# minimap2-sys 0.1.30 with minimap2 2.30-r1287

The crate minimap2-sys 0.1.30+minimap2.2.30 on crates.io carries minimap2 2.29-r1283. Its C sources here are replaced
by those of the minimap2 v2.30 release (r1287), which scales minisplice splice-site scores on loading (`--spsc-scale`,
0.7 by default; `mm_idx_spsc_read2`). The public structs are unchanged between the two, so the crate's generated
bindings still describe them; IsoForge declares `mm_idx_spsc_read2` itself. With these sources, IsoForge's in-process
alignment gives the records of `minimap2` 2.30-r1287 exactly. Used through `[patch.crates-io]` in IsoForge's Cargo.toml.
