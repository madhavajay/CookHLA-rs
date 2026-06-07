//! CONVERT_IN (target side) — turn the QC'd target (`MHC.QC` PED/DAT/BIM) into the `MHC.QC.vcf`
//! that IMPUTE feeds to Beagle. Native port of `HLA_Imputation_BEAGLE5.CONVERT_IN`'s target path:
//!
//! PED+DAT → bgl (`linkage2beagle`) → refine marker positions to the reference
//! (`excluding_snp_and_refine_target_position-*.R`, **including its first-matched-marker
//! off-by-one bug**, which we preserve) → subset (`Panel_subset.py`) → GC trick → `beagle2vcf`.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::Result;

use super::qc::QcOutputs;
use crate::bgl::{Bgl, BglLine, MarkerRecord, Markers};
use crate::convert::{beagle2vcf, linkage2beagle, parse_dat};

/// Port of `excluding_snp_and_refine_target_position-v1COOK02222017.R`.
///
/// For each target marker (in file order), if it exists in the reference, take the **reference**
/// row (id, ref-bp, ref alleles). Then — preserving the R script's off-by-one bug — drop the
/// *first* matched marker, and sort the rest by reference base position.
pub fn refine_target_markers(target: &Markers, reference: &Markers) -> Markers {
    let ref_by_id: HashMap<&str, &MarkerRecord> = reference
        .records
        .iter()
        .map(|r| (r.id.as_str(), r))
        .collect();

    let mut matched: Vec<MarkerRecord> = target
        .records
        .iter()
        .filter_map(|t| ref_by_id.get(t.id.as_str()).map(|r| (*r).clone()))
        .collect();

    // R bug: the first matched marker is written to index 0 (a no-op in 1-based R) and lost.
    if !matched.is_empty() {
        matched.remove(0);
    }
    // R: order(new_data_marker[,2]) — stable sort by base position.
    matched.sort_by_key(|r| r.bp);

    Markers { records: matched }
}

/// Port of `Panel_subset.py` for the `(indv = all, markers = <set>)` case used by CONVERT_IN:
/// keep only the bgl `M` rows / marker rows whose id is in `selected`. Returns a bgl whose `M`
/// rows are reordered to match `refined_markers` so `beagle2vcf` pairs them correctly.
fn panel_subset(bgl: &Bgl, refined_markers: &Markers, selected: &HashSet<&str>) -> Bgl {
    // Index the bgl's marker rows by id.
    let by_id: HashMap<&str, &Vec<String>> = bgl
        .lines
        .iter()
        .filter_map(|l| match l {
            BglLine::Marker { label, alleles } if selected.contains(label.as_str()) => {
                Some((label.as_str(), alleles))
            }
            _ => None,
        })
        .collect();

    let mut lines: Vec<BglLine> = bgl
        .lines
        .iter()
        .filter(|l| matches!(l, BglLine::Header(_)))
        .cloned()
        .collect();
    for m in &refined_markers.records {
        if let Some(alleles) = by_id.get(m.id.as_str()) {
            lines.push(BglLine::Marker {
                label: m.id.clone(),
                alleles: (*alleles).clone(),
            });
        }
    }
    Bgl { lines }
}

/// Run CONVERT_IN (target side): produce the imputation target VCF text (`MHC.QC.vcf`).
///
/// `reference_markers_path` is the panel `.markers`; `chrom` is the VCF chromosome label (`6`).
pub fn convert_in_target(
    qc: &QcOutputs,
    reference_markers_path: &Path,
    chrom: &str,
) -> Result<String> {
    // 1. PED + DAT → bgl.
    let ped = std::fs::read_to_string(&qc.nopheno_ped)?;
    let dat = std::fs::read_to_string(&qc.dat)?;
    let marker_ids = parse_dat(&dat);
    let bgl_text = linkage2beagle(&ped, &marker_ids)?;
    let bgl = Bgl::parse(&bgl_text)?;

    // 2. Refine target marker positions against the reference (with the preserved R bug).
    let target_markers = markers_from_bim(&qc.bim)?;
    let reference_markers = Markers::read(reference_markers_path)?;
    let refined_markers = refine_target_markers(&target_markers, &reference_markers);

    // 3. Subset the bgl to the refined marker set, aligned to the refined marker order.
    let selected: HashSet<&str> = refined_markers
        .records
        .iter()
        .map(|r| r.id.as_str())
        .collect();
    let refined_bgl = panel_subset(&bgl, &refined_markers, &selected);

    // 4. GC trick, then 5. beagle2vcf.
    let (gc_bgl, gc_markers) = refined_bgl.gc_trick(&refined_markers)?;
    beagle2vcf(chrom, &gc_markers, &gc_bgl, "0")
}

/// Build a [`Markers`] table from a `.bim`'s (id, bp, a1, a2) columns, in file (target) order.
fn markers_from_bim(bim_path: &Path) -> Result<Markers> {
    let bim = crate::plink::Bim::read(bim_path)?;
    Ok(Markers {
        records: bim
            .records
            .iter()
            .map(|r| MarkerRecord {
                id: r.id.clone(),
                bp: r.bp,
                a1: r.a1.clone(),
                a2: r.a2.clone(),
            })
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refine_preserves_first_match_dropping_bug() {
        // target order: t1 (in ref), t2 (in ref), t3 (not in ref). Ref has t1@200, t2@100.
        let target = Markers::parse("t1 999 A G\nt2 999 C T\nt3 999 A C\n").unwrap();
        let reference = Markers::parse("t1 200 A G\nt2 100 C T\n").unwrap();
        let refined = refine_target_markers(&target, &reference);
        // t1 (first match) dropped by the bug; only t2 remains (ref position 100).
        assert_eq!(refined.records.len(), 1);
        assert_eq!(refined.records[0].id, "t2");
        assert_eq!(refined.records[0].bp, 100);
    }

    #[test]
    fn refine_sorts_remaining_by_ref_bp() {
        // first match dropped, the rest sorted by reference bp.
        let target = Markers::parse("d 9 A G\na 9 A G\nb 9 A G\nc 9 A G\n").unwrap();
        let reference = Markers::parse("a 400 A G\nb 100 A G\nc 300 A G\nd 500 A G\n").unwrap();
        let refined = refine_target_markers(&target, &reference);
        // matched in target order [d,a,b,c]; drop d; remaining [a,b,c] sorted by bp → b(100),c(300),a(400).
        let ids: Vec<&str> = refined.records.iter().map(|r| r.id.as_str()).collect();
        assert_eq!(ids, ["b", "c", "a"]);
    }
}
