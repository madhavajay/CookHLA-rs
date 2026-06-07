//! Beagle v3-style companion formats CookHLA shuttles between PLINK, MACH, and Beagle:
//! the `.markers` table and the `.bgl(.phased)` genotype file, plus the "GC trick".
//!
//! The GC trick (`src/bgl2GC_trick_bgl.py`) temporarily recodes every marker's two alleles to
//! `G`/`C` because Beagle's utilities choke on CookHLA's `P`/`A` (presence/absence) allele
//! names. It is applied to both target and reference before conversion to VCF, and reversed on
//! the way out. Porting it natively removes a Python hop from the hot path.

pub mod markers;
pub mod phased;

pub use markers::{MarkerRecord, Markers};
pub use phased::{Bgl, BglLine};
