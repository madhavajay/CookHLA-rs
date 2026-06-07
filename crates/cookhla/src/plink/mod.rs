//! PLINK file formats used by CookHLA: `.bim`, `.fam`, and `.frq` (from `plink --freq`).
//!
//! CookHLA invokes the `plink` binary for the heavy operations (`--make-bed`, `--recode`,
//! `--freq`, `--flip`, ...). We keep that binary for now (it is fast and compiled), but we
//! must read and write these text companion files ourselves to do the QC/merge logic that
//! CookHLA does in awk/sed/Perl/pandas. The `.bed` genotype matrix is handled separately
//! (Phase 2 continuation) because most CookHLA QC works off `.bim`/`.fam`/`.frq`.

pub mod bim;
pub mod fam;
pub mod frq;

pub use bim::{Bim, BimRecord};
pub use fam::{Fam, FamRecord};
pub use frq::{Frq, FrqRecord};
