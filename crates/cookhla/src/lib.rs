//! CookHLA-rs — a Rust port of [CookHLA](https://github.com/WansonChoi/CookHLA),
//! an accurate and efficient HLA imputation method.
//!
//! The original Python/R/csh/Perl pipeline lives in the `repos/CookHLA` submodule and is the
//! reference implementation. This crate ports it bottom-up: file formats first, then the
//! converters and QC, the consensus caller, the local-embedding panels, and finally the
//! imputation driver (which calls `beagle-rs` in-process). See `TODO.md` for the phased plan.
//!
//! Parity goal (v1): identical HLA *calls* to the Python reference, probabilities within a
//! small epsilon. Speed goal: replace the interpreted orchestration/consensus with compiled
//! Rust while keeping the already-fast PLINK/MACH binaries (for now).

pub mod bgl;
pub mod consensus;
pub mod convert;
pub mod front;
pub mod impute;
pub mod pipeline;
pub mod plink;

/// The eight classical HLA genes CookHLA imputes, in the order the reference uses.
pub const HLA_NAMES: [&str; 8] = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"];

/// Split a line on runs of ASCII whitespace, skipping empties (mirrors `split /\s+/` / awk).
pub(crate) fn ws_fields(line: &str) -> impl Iterator<Item = &str> {
    line.split_ascii_whitespace()
}
