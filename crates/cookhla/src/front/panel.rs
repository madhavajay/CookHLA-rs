//! Reference panel building (`Make_EXON234_Panel` + `HLA_MultipleRefs.Make_ExonN_Panel`) — the
//! local-embedding trick. From the reference panel, build the per-exon phased VCFs that IMPUTE
//! uses as `ref=`.
//!
//! Steps (native ports of the Python/R glue):
//! 1. Keep SNPs + 4-digit HLA markers (drop `AA_`/`SNP_`/`INS_`/2-digit HLA).
//! 2. `HLA2EXON234`: triple each HLA marker into `_exon2/_exon3/_exon4` at the exon midpoint
//!    positions (class II — DRB1/DQ/DP — has no exon 4).
//! 3. `redefineBP` (disperse duplicate positions) + sort by position; reorder the bgl to match.
//! 4. Per exon: subset to that exon's HLA markers + SNPs → GC trick → `beagle2vcf` → phase (`/`→`|`).

use std::collections::{HashMap, HashSet};

use anyhow::Result;

use crate::bgl::{Bgl, BglLine, MarkerRecord, Markers};
use crate::convert::beagle2vcf;

/// HLA exon-midpoint base positions (hg18), from `Make_EXON234_Panel.HLA2EXON234`.
const EXON2_POS: &[(&str, &str)] = &[
    ("A", "30018647"),
    ("C", "31347489"),
    ("B", "31432578"),
    ("DRB1", "32659998"),
    ("DQA1", "32717189"),
    ("DQB1", "32740687"),
    ("DPA1", "33145518"),
    ("DPB1", "33156558"),
];
const EXON3_POS: &[(&str, &str)] = &[
    ("A", "30019161"),
    ("C", "31346966"),
    ("B", "31432060"),
    ("DRB1", "32657452"),
    ("DQA1", "32717867"),
    ("DQB1", "32737862"),
    ("DPA1", "33144914"),
    ("DPB1", "33160845"),
];
const EXON4_POS: &[(&str, &str)] = &[("A", "30020015"), ("C", "31346103"), ("B", "31431210")];

/// The gene token of an HLA marker label: `HLA_DRB1_1501` → `DRB1`.
fn gene_of(label: &str) -> Option<&str> {
    label.strip_prefix("HLA_").and_then(|r| r.split('_').next())
}

fn pos_lookup(table: &[(&str, &str)], gene: &str) -> Option<i64> {
    table
        .iter()
        .find(|(g, _)| *g == gene)
        .and_then(|(_, p)| p.parse().ok())
}

/// STEP1: keep SNPs and 4-digit HLA markers; drop `AA_`/`SNP_`/`INS_` and 2-digit HLA
/// (`HLA_<gene>_<2 digits>`).
fn select_markers(ref_markers: &Markers) -> HashSet<String> {
    ref_markers
        .records
        .iter()
        .filter(|m| !excluded(&m.id))
        .map(|m| m.id.clone())
        .collect()
}

fn excluded(id: &str) -> bool {
    if id.starts_with("AA_") || id.starts_with("SNP_") || id.starts_with("INS_") {
        return true;
    }
    if id.starts_with("HLA_") {
        // 2-digit HLA: last `_`-segment is exactly two ASCII digits.
        if let Some(last) = id.rsplit('_').next() {
            if last.len() == 2 && last.bytes().all(|b| b.is_ascii_digit()) {
                return true;
            }
        }
    }
    false
}

/// `HLA2EXON234`: triple HLA markers into per-exon copies; SNPs pass through once.
fn hla2exon234(select: &HashSet<String>, ref_bgl: &Bgl, ref_markers: &Markers) -> (Bgl, Markers) {
    // --- bgl ---
    let mut lines = Vec::new();
    for l in &ref_bgl.lines {
        match l {
            BglLine::Header(_) => lines.push(l.clone()),
            BglLine::Marker { label, alleles } => {
                if !select.contains(label) {
                    continue;
                }
                if label.starts_with("HLA_") {
                    let class_i = matches!(gene_of(label), Some("A" | "B" | "C"));
                    for (exon, on) in [(2u8, true), (3, true), (4, class_i)] {
                        if on {
                            lines.push(BglLine::Marker {
                                label: format!("{label}_exon{exon}"),
                                alleles: alleles.clone(),
                            });
                        }
                    }
                } else {
                    lines.push(l.clone());
                }
            }
        }
    }

    // --- markers ---
    let mut recs = Vec::new();
    for m in &ref_markers.records {
        if !select.contains(&m.id) {
            continue;
        }
        if m.id.starts_with("HLA_") {
            let gene = gene_of(&m.id).unwrap_or("");
            for (exon, table) in [(2u8, EXON2_POS), (3, EXON3_POS), (4, EXON4_POS)] {
                if let Some(bp) = pos_lookup(table, gene) {
                    recs.push(MarkerRecord {
                        id: format!("{}_exon{exon}", m.id),
                        bp,
                        a1: m.a1.clone(),
                        a2: m.a2.clone(),
                    });
                }
            }
        } else {
            recs.push(m.clone());
        }
    }

    (Bgl { lines }, Markers { records: recs })
}

/// `redefineBP`: push each duplicate base position up by 1 until free, then sort by position.
fn redefine_bp(markers: &Markers) -> Markers {
    let mut occupied: HashSet<i64> = HashSet::new();
    let mut placed: Vec<MarkerRecord> = Vec::with_capacity(markers.records.len());
    for m in &markers.records {
        let mut bp = m.bp;
        while occupied.contains(&bp) {
            bp += 1;
        }
        occupied.insert(bp);
        placed.push(MarkerRecord { bp, ..m.clone() });
    }
    placed.sort_by_key(|m| m.bp);
    Markers { records: placed }
}

/// `BGL2SortBGL_WS`: reorder the bgl `M` rows to the marker order (headers kept in place).
fn bgl_sort(order: &Markers, bgl: &Bgl) -> Bgl {
    let by_id: HashMap<&str, &Vec<String>> = bgl
        .lines
        .iter()
        .filter_map(|l| match l {
            BglLine::Marker { label, alleles } => Some((label.as_str(), alleles)),
            _ => None,
        })
        .collect();
    let mut lines: Vec<BglLine> = bgl
        .lines
        .iter()
        .filter(|l| matches!(l, BglLine::Header(_)))
        .cloned()
        .collect();
    for m in &order.records {
        if let Some(alleles) = by_id.get(m.id.as_str()) {
            lines.push(BglLine::Marker {
                label: m.id.clone(),
                alleles: (*alleles).clone(),
            });
        }
    }
    Bgl { lines }
}

/// Build the exon234 panel (bgl + markers) from a reference panel.
pub fn make_exon234(ref_bgl: &Bgl, ref_markers: &Markers) -> (Bgl, Markers) {
    let select = select_markers(ref_markers);
    let (e_bgl, e_markers) = hla2exon234(&select, ref_bgl, ref_markers);
    let sorted_markers = redefine_bp(&e_markers);
    let sorted_bgl = bgl_sort(&sorted_markers, &e_bgl);
    (sorted_bgl, sorted_markers)
}

/// Build one exon's phased reference VCF (`HM_CEU_REF.exon<N>.phased.vcf`) from the exon234 panel.
pub fn make_exon_panel_vcf(
    exon: u8,
    exon234_bgl: &Bgl,
    exon234_markers: &Markers,
) -> Result<String> {
    let suffix = format!("_exon{exon}");
    let keep = |label: &str| !label.starts_with("HLA_") || label.ends_with(&suffix);

    // Subset markers + bgl to this exon's HLA markers + the SNPs.
    let markers = Markers {
        records: exon234_markers
            .records
            .iter()
            .filter(|m| keep(&m.id))
            .cloned()
            .collect(),
    };
    let kept_ids: HashSet<&str> = markers.records.iter().map(|m| m.id.as_str()).collect();
    let bgl = Bgl {
        lines: exon234_bgl
            .lines
            .iter()
            .filter(|l| match l {
                BglLine::Header(_) => true,
                BglLine::Marker { label, .. } => kept_ids.contains(label.as_str()),
            })
            .cloned()
            .collect(),
    };

    // GC trick → beagle2vcf → phase (`/` → `|`, since the reference is phased).
    let (gc_bgl, gc_markers) = bgl.gc_trick(&markers)?;
    let vcf = beagle2vcf("6", &gc_markers, &gc_bgl, "0")?;
    Ok(vcf.replace('/', "|"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn excludes_aa_snp_ins_and_2digit_hla() {
        assert!(excluded("AA_A_-15_30018338"));
        assert!(excluded("SNP_A_30018337"));
        assert!(excluded("HLA_A_01")); // 2-digit
        assert!(!excluded("HLA_A_0101")); // 4-digit kept
        assert!(!excluded("rs7772982")); // SNP kept
    }

    #[test]
    fn gene_extraction() {
        assert_eq!(gene_of("HLA_DRB1_1501"), Some("DRB1"));
        assert_eq!(gene_of("HLA_A_0101"), Some("A"));
    }

    #[test]
    fn redefine_bp_disperses_duplicates() {
        let m = Markers::parse("a 100 P A\nb 100 P A\nc 100 P A\n").unwrap();
        let r = redefine_bp(&m);
        let bps: Vec<i64> = r.records.iter().map(|x| x.bp).collect();
        assert_eq!(bps, [100, 101, 102]); // dispersed then sorted
    }

    #[test]
    fn class_ii_has_no_exon4() {
        let select: HashSet<String> = ["HLA_DRB1_1501".to_string()].into_iter().collect();
        let bgl = Bgl::parse("I id s s\nM HLA_DRB1_1501 P A\n").unwrap();
        let markers = Markers::parse("HLA_DRB1_1501 32000000 P A\n").unwrap();
        let (e_bgl, e_markers) = hla2exon234(&select, &bgl, &markers);
        let labels: Vec<&str> = e_bgl.marker_labels().collect();
        assert_eq!(labels, ["HLA_DRB1_1501_exon2", "HLA_DRB1_1501_exon3"]); // no exon4
        assert_eq!(e_markers.records.len(), 2);
    }
}
