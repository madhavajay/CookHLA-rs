//! Parity for reference panel building against the golden per-exon phased VCFs. Builds the
//! exon234 panel from the example reference and checks each exon's phased VCF (per-marker
//! REF/ALT/GT) against the golden `HM_CEU_REF.exon<N>.phased.vcf`.

use std::collections::HashMap;
use std::path::PathBuf;

use cookhla::bgl::{Bgl, Markers};
use cookhla::front::panel::{make_exon234, make_exon_panel_vcf};

fn manifest() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}
fn golden_dir() -> PathBuf {
    manifest()
        .join("../../fixtures/example/golden")
        .canonicalize()
        .unwrap_or_else(|_| PathBuf::from("/nonexistent"))
}

fn parse_vcf(text: &str) -> (Vec<String>, HashMap<String, String>) {
    let mut samples = Vec::new();
    let mut rows = HashMap::new();
    for line in text.lines() {
        if line.starts_with("##") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if line.starts_with("#CHROM") {
            samples = f[9..].iter().map(|s| s.to_string()).collect();
            continue;
        }
        if f.len() < 10 {
            continue;
        }
        rows.insert(
            f[2].to_string(),
            format!("{} {} {}", f[3], f[4], f[9..].join(" ")),
        );
    }
    (samples, rows)
}

#[test]
fn panel_building_matches_golden() {
    let g = golden_dir();
    let ref_prefix = manifest().join("../../repos/CookHLA/example/HM_CEU_REF");
    let golden_e2 = g.join("HM_CEU_REF.exon2.phased.vcf");
    if !golden_e2.exists() {
        eprintln!(
            "SKIP: golden panels not at {}. Run docker/oracle-run.sh.",
            g.display()
        );
        return;
    }

    let ref_bgl = Bgl::read(ref_prefix.with_extension("bgl.phased")).unwrap();
    let ref_markers = Markers::read(ref_prefix.with_extension("markers")).unwrap();
    let (e234_bgl, e234_markers) = make_exon234(&ref_bgl, &ref_markers);
    eprintln!("exon234 panel: {} markers", e234_markers.records.len());

    for exon in 2u8..=4 {
        let golden_path = g.join(format!("HM_CEU_REF.exon{exon}.phased.vcf"));
        let got = make_exon_panel_vcf(exon, &e234_bgl, &e234_markers).unwrap();
        let (got_s, got_rows) = parse_vcf(&got);
        let (want_s, want_rows) = parse_vcf(&std::fs::read_to_string(&golden_path).unwrap());
        assert_eq!(got_s, want_s, "exon{exon}: sample order differs");
        eprintln!(
            "exon{exon}: rust {} rows, golden {} rows",
            got_rows.len(),
            want_rows.len()
        );
        assert_eq!(
            got_rows.len(),
            want_rows.len(),
            "exon{exon}: marker count differs"
        );
        let mut mism = 0;
        for (id, want) in &want_rows {
            if got_rows.get(id) != Some(want) {
                if mism < 8 {
                    eprintln!(
                        "exon{exon} row mismatch {id}: rust {:?} vs golden {:?}",
                        got_rows.get(id),
                        want
                    );
                }
                mism += 1;
            }
        }
        assert_eq!(mism, 0, "exon{exon}: phased VCF rows must match golden");
        eprintln!(
            "exon{exon} panel parity OK: {} markers match golden",
            want_rows.len()
        );
    }
}
