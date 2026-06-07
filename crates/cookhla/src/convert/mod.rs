//! Native replacements for the small Beagle utility JARs CookHLA shells out to.
//!
//! These are pure format converters; porting them removes five `java -jar` hops from the
//! pipeline. Each is verified against golden output captured from the original jar
//! (`repos/CookHLA/dependency/*.jar`).
//!
//! - [`beagle2vcf`] — `.bgl` (+ `.markers`) → VCF (replaces `beagle2vcf.jar`).
//! - [`linkage2beagle`] — `.ped` + `.dat` → `.bgl` (replaces `linkage2beagle.jar`).
//! - [`beagle2linkage`] — `.bgl` → `.ped` + `.dat` (replaces `beagle2linkage.jar`).
//!
//! Still to port (Phase 3): `vcf2beagle`, `transpose`.

pub mod beagle2linkage;
pub mod beagle2vcf;
pub mod linkage2beagle;

pub use beagle2linkage::{beagle2linkage, LinkageOut};
pub use beagle2vcf::beagle2vcf;
pub use linkage2beagle::{linkage2beagle, parse_dat};
