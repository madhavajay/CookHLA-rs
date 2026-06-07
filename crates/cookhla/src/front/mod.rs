//! The "front half" of CookHLA: turn raw target genotypes + a reference panel into the
//! CONVERT_IN inputs the (already-ported, verified) IMPUTE→CONSENSUS back half consumes.
//!
//! These stages orchestrate the fast `plink` binary (kept, per the strategy) plus native Rust
//! replacements for the awk/sed/`merge_tables.pl`/grep glue:
//! - [`plink`] — a thin wrapper over the `plink` executable.
//! - [`tables`] — `merge_tables.pl` (index join) + the small table ops the QC uses.
//! - [`qc`] — `EXTRACT_MHC` + `FLIP` (strand flip, frequency QC, recode) → `MHC.QC.*`.

pub mod agm;
pub mod convert_in;
pub mod fixinput;
pub mod liftover;
pub mod panel;
pub mod plink;
pub mod qc;
pub mod tables;

pub use plink::Plink;
