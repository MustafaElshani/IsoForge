# rust-htslib 0.47.1 with the bindings shipped by hts-sys

The crate rust-htslib 0.47.1 on crates.io generates its htslib bindings at build time (`hts-sys` with the `bindgen`
feature), which needs libclang on the build machine. Here that feature is left out, so the bindings shipped with
`hts-sys` are used and the build needs a C compiler only. The integer conversions and one field name (`isize`) in
`src/bam`, `src/bcf`, `src/bgzf` and `src/faidx` follow the shipped bindings. Nothing else differs from the published
crate. Used through `[patch.crates-io]` in IsoForge's Cargo.toml.
