//! Adaptive genetic map (AGM) exon-splitting — `Make_EXON234_AGM` + `HLA_MultipleRefs.Make_ExonN_AGM`.
//! Turn the provided whole-region map (`*.mach_step.avg.clpsB`, SNP genetic distances) into the
//! per-exon maps IMPUTE passes to Beagle as `map=`.
//!
//! The HLA exon markers aren't in the SNP map, so they're placed at the genetic midpoint of their
//! two flanking SNPs (then nudged apart by ε so positions stay strictly increasing) — the
//! `GEN_stitch_GD` step. Beagle only needs a monotone map close to the original, so small
//! floating-point rendering differences are immaterial to the final calls.

use std::collections::HashMap;

use anyhow::{Context, Result};
use std::path::Path;

use crate::bgl::Markers;

const EPS: f64 = 1e-12;

/// One genetic-map row: `chrom  id  cM  bp` (Beagle/PLINK map format).
#[derive(Debug, Clone)]
pub struct AgmRow {
    pub chr: String,
    pub id: String,
    pub gd: f64,
    pub bp: i64,
}

/// Parse a `.clpsB` / AGM file (`chr id gd bp`, whitespace-separated).
pub fn parse_agm(text: &str) -> Result<Vec<AgmRow>> {
    let mut rows = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 4 {
            continue;
        }
        rows.push(AgmRow {
            chr: f[0].to_string(),
            id: f[1].to_string(),
            gd: f[2]
                .parse()
                .with_context(|| format!("agm: bad GD on line {}", i + 1))?,
            bp: f[3]
                .parse()
                .with_context(|| format!("agm: bad bp on line {}", i + 1))?,
        });
    }
    Ok(rows)
}

/// Serialize AGM rows back to `chr\tid\tgd\tbp` text.
pub fn agm_to_text(rows: &[AgmRow]) -> String {
    let mut s = String::new();
    for r in rows {
        s.push_str(&format!("{}\t{}\t{}\t{}\n", r.chr, r.id, r.gd, r.bp));
    }
    s
}

fn is_snp_marker(id: &str) -> bool {
    !(id.starts_with("AA_")
        || id.starts_with("SNP_")
        || id.starts_with("HLA_")
        || id.starts_with("INS_"))
}

/// `Make_EXON234_AGM`: merge the SNP genetic map with the exon234 panel markers and interpolate
/// the HLA exon markers' genetic distances.
pub fn make_exon234_agm(provided: &[AgmRow], exon234_markers: &Markers) -> Vec<AgmRow> {
    // SNP genetic distances from the provided map.
    let gd_of: HashMap<&str, f64> = provided
        .iter()
        .filter(|r| is_snp_marker(&r.id))
        .map(|r| (r.id.as_str(), r.gd))
        .collect();
    let provided_bp: HashMap<&str, i64> = provided.iter().map(|r| (r.id.as_str(), r.bp)).collect();

    // Outer union of (provided SNPs) and (exon234 markers), keyed by id.
    let panel_bp: HashMap<&str, i64> = exon234_markers
        .records
        .iter()
        .map(|m| (m.id.as_str(), m.bp))
        .collect();

    let mut seen = std::collections::HashSet::new();
    let mut rows: Vec<AgmRow> = Vec::new();
    let push = |id: &str, rows: &mut Vec<AgmRow>, seen: &mut std::collections::HashSet<String>| {
        if !seen.insert(id.to_string()) {
            return;
        }
        let bp = panel_bp
            .get(id)
            .copied()
            .or_else(|| provided_bp.get(id).copied())
            .unwrap_or(0);
        let gd = gd_of.get(id).copied().unwrap_or(0.0);
        rows.push(AgmRow {
            chr: "6".into(),
            id: id.to_string(),
            gd,
            bp,
        });
    };
    // Panel markers first, then any provided SNPs not in the panel.
    for m in &exon234_markers.records {
        push(&m.id, &mut rows, &mut seen);
    }
    for r in provided.iter().filter(|r| is_snp_marker(&r.id)) {
        push(&r.id, &mut rows, &mut seen);
    }

    rows.sort_by_key(|r| r.bp);
    stitch_gd(&mut rows);
    rows
}

/// `GEN_stitch_GD`: assign each run of HLA markers (genetic distance 0) the midpoint of the two
/// flanking SNP distances, nudged by ε per marker so the map stays strictly increasing.
fn stitch_gd(rows: &mut [AgmRow]) {
    let n = rows.len();
    if n == 0 {
        return;
    }
    let gd: Vec<f64> = rows.iter().map(|r| r.gd).collect();
    let mut new_gd: Vec<f64> = Vec::with_capacity(n);
    new_gd.push(0.0); // first row, per the R code

    let mut i = 1;
    while i < n {
        if gd[i - 1] != 0.0 && gd[i] == 0.0 {
            let start = i;
            let start_cap = start - 1;
            let mut end = i + 1;
            while end < n && gd[end] == 0.0 {
                end += 1;
            }
            // `end` now points at the first non-zero after the chunk (the trailing flank).
            let end_cap = end.min(n - 1);
            let last_hla = end - 1;

            let mid = (gd[start_cap] + gd[end_cap]) / 2.0;
            let mut acc = mid;
            new_gd.push(mid); // first HLA marker in the chunk
            let mut j = start + 1;
            while j <= last_hla {
                acc += EPS;
                new_gd.push(acc);
                j += 1;
            }
            i = last_hla; // R sets i = end (its `end` is the last HLA index)
        } else {
            new_gd.push(gd[i]);
        }
        i += 1;
    }

    for (r, g) in rows.iter_mut().zip(new_gd) {
        r.gd = g;
    }
}

/// `Make_ExonN_AGM`: subset the exon234 AGM to SNPs + this exon's HLA markers (`*_exon<N>`).
pub fn make_exon_agm(exon: u8, exon234_agm: &[AgmRow]) -> Vec<AgmRow> {
    let suffix = format!("_exon{exon}");
    exon234_agm
        .iter()
        .filter(|r| !r.id.starts_with("HLA_") || r.id.ends_with(&suffix))
        .cloned()
        .collect()
}

/// Convenience: read a provided AGM file and build the exon234 + per-exon maps, writing the
/// per-exon maps to `out_dir/<stem>.exon<N>.txt`.
pub fn build_exon_maps(
    provided_path: &Path,
    exon234_markers: &Markers,
) -> Result<HashMap<u8, Vec<AgmRow>>> {
    let provided = parse_agm(&std::fs::read_to_string(provided_path)?)?;
    let e234 = make_exon234_agm(&provided, exon234_markers);
    Ok([2u8, 3, 4]
        .into_iter()
        .map(|e| (e, make_exon_agm(e, &e234)))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snp_classification() {
        assert!(is_snp_marker("rs123"));
        assert!(!is_snp_marker("HLA_A_0101_exon2"));
        assert!(!is_snp_marker("AA_A_-15"));
    }

    #[test]
    fn stitch_interpolates_hla_run_to_midpoint() {
        // SNP(gd 1.0) | HLA HLA | SNP(gd 2.0): HLA markers get midpoint 1.5, 1.5+eps.
        let mut rows = vec![
            AgmRow {
                chr: "6".into(),
                id: "rsL".into(),
                gd: 1.0,
                bp: 10,
            },
            AgmRow {
                chr: "6".into(),
                id: "HLA_A_0101_exon2".into(),
                gd: 0.0,
                bp: 20,
            },
            AgmRow {
                chr: "6".into(),
                id: "HLA_A_0201_exon2".into(),
                gd: 0.0,
                bp: 21,
            },
            AgmRow {
                chr: "6".into(),
                id: "rsR".into(),
                gd: 2.0,
                bp: 30,
            },
        ];
        stitch_gd(&mut rows);
        assert_eq!(rows[0].gd, 0.0); // first row forced to 0 by the R code
        assert!((rows[1].gd - 1.5).abs() < 1e-15);
        assert!((rows[2].gd - (1.5 + EPS)).abs() < 1e-15);
        assert!((rows[3].gd - 2.0).abs() < 1e-15);
    }

    #[test]
    fn exon_subset_keeps_snps_and_only_that_exon() {
        let e234 = vec![
            AgmRow {
                chr: "6".into(),
                id: "rs1".into(),
                gd: 0.1,
                bp: 1,
            },
            AgmRow {
                chr: "6".into(),
                id: "HLA_A_0101_exon2".into(),
                gd: 0.2,
                bp: 2,
            },
            AgmRow {
                chr: "6".into(),
                id: "HLA_A_0101_exon3".into(),
                gd: 0.3,
                bp: 3,
            },
        ];
        let e2 = make_exon_agm(2, &e234);
        let ids: Vec<&str> = e2.iter().map(|r| r.id.as_str()).collect();
        assert_eq!(ids, ["rs1", "HLA_A_0101_exon2"]);
    }
}
