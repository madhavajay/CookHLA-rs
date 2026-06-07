//! Parity for the QC stage against the golden `MHC.QC.*`. Runs `run_qc` on the golden
//! (de-ambiguated) target + the example reference and checks the marker selection, the
//! reference-repositioned `.bim`, and the `.dat`. Skips loudly if plink or the golden are absent.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use cookhla::front::{convert_in::convert_in_target, qc::run_qc, Plink};
use cookhla::plink::Bim;

fn golden_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/example/golden")
        .canonicalize()
        .unwrap_or_else(|_| PathBuf::from("/nonexistent"))
}

fn example_ref() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../repos/CookHLA/example/HM_CEU_REF")
}

#[test]
fn qc_matches_golden_marker_selection() {
    let g = golden_dir();
    let noambig = g.join("1958BC.hg19.COPY.LiftDown_hg18.NoAmbig.bed");
    let golden_dat = g.join("1958BC+HM_CEU_REF.MHC.QC.dat");
    let golden_bim = g.join("1958BC+HM_CEU_REF.MHC.QC.bim");
    let Some(plink) = Plink::locate() else {
        eprintln!("SKIP: plink not found (set $PLINK or vendor repos/CookHLA/dependency/plink)");
        return;
    };
    if !noambig.exists() || !golden_dat.exists() {
        eprintln!(
            "SKIP: golden QC inputs not at {}. Run docker/oracle-run.sh.",
            g.display()
        );
        return;
    }

    let tmp = std::env::temp_dir().join("cookhla_qc_parity");
    std::fs::create_dir_all(&tmp).unwrap();
    let input_prefix = g.join("1958BC.hg19.COPY.LiftDown_hg18.NoAmbig");
    let mhc_prefix = tmp.join("1958BC+HM_CEU_REF.MHC");

    let out = run_qc(
        &plink,
        &input_prefix,
        &example_ref(),
        &mhc_prefix,
        /*small_sample=*/ true,
    )
    .expect("run_qc");

    // --- marker selection: the .dat marker list (set + order) must match golden ---
    let got_dat = std::fs::read_to_string(&out.dat).unwrap();
    let want_dat = std::fs::read_to_string(&golden_dat).unwrap();
    let got_markers: Vec<&str> = got_dat.lines().collect();
    let want_markers: Vec<&str> = want_dat.lines().collect();
    eprintln!(
        "QC markers: rust {}, golden {}",
        got_markers.len(),
        want_markers.len()
    );
    assert_eq!(
        got_markers.len(),
        want_markers.len(),
        "marker count differs"
    );
    let got_set: std::collections::BTreeSet<_> = got_markers.iter().collect();
    let want_set: std::collections::BTreeSet<_> = want_markers.iter().collect();
    let missing: Vec<_> = want_set.difference(&got_set).take(10).collect();
    let extra: Vec<_> = got_set.difference(&want_set).take(10).collect();
    assert!(
        missing.is_empty(),
        "markers missing vs golden (first 10): {missing:?}"
    );
    assert!(
        extra.is_empty(),
        "extra markers vs golden (first 10): {extra:?}"
    );

    // --- repositioned .bim: every selected marker must land at the golden's reference position ---
    let got_bim = Bim::read(&out.bim).unwrap();
    let want_bim = Bim::read(Path::new(&golden_bim)).unwrap();
    let want_pos: std::collections::HashMap<&str, i64> = want_bim
        .records
        .iter()
        .map(|r| (r.id.as_str(), r.bp))
        .collect();
    let mut pos_mismatch = 0;
    for r in &got_bim.records {
        if let Some(&wp) = want_pos.get(r.id.as_str()) {
            if wp != r.bp {
                if pos_mismatch < 10 {
                    eprintln!("POS mismatch {}: rust {} vs golden {}", r.id, r.bp, wp);
                }
                pos_mismatch += 1;
            }
        }
    }
    assert_eq!(pos_mismatch, 0, "reference repositioning mismatches");
    eprintln!(
        "QC parity OK: {} markers, positions match golden",
        got_markers.len()
    );

    // --- CONVERT_IN: produce MHC.QC.vcf and compare to golden (per-marker REF/ALT/GT) ---
    let ref_markers = example_ref().with_extension("markers");
    let got_vcf = convert_in_target(&out, &ref_markers, "6").expect("convert_in_target");
    let golden_vcf = std::fs::read_to_string(g.join("1958BC+HM_CEU_REF.MHC.QC.vcf")).unwrap();

    let (got_samples, got_rows) = parse_vcf(&got_vcf);
    let (want_samples, want_rows) = parse_vcf(&golden_vcf);
    assert_eq!(got_samples, want_samples, "MHC.QC.vcf sample order differs");
    eprintln!(
        "CONVERT_IN markers: rust {}, golden {}",
        got_rows.len(),
        want_rows.len()
    );
    assert_eq!(
        got_rows.len(),
        want_rows.len(),
        "MHC.QC.vcf marker count differs"
    );
    let mut row_mismatch = 0;
    for (id, want) in &want_rows {
        match got_rows.get(id) {
            Some(got) if got == want => {}
            other => {
                if row_mismatch < 10 {
                    eprintln!(
                        "VCF row mismatch {id}: rust {:?} vs golden {:?}",
                        other, want
                    );
                }
                row_mismatch += 1;
            }
        }
    }
    assert_eq!(row_mismatch, 0, "MHC.QC.vcf rows must match golden");
    eprintln!(
        "CONVERT_IN parity OK: {} marker rows match golden MHC.QC.vcf",
        want_rows.len()
    );
}

/// Parse a beagle2vcf VCF into (sample names, {marker_id -> "REF ALT gt1 gt2 ..."}).
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
        // key on REF, ALT, and the genotype columns (skip POS so a benign coordinate tie can't fail it).
        let payload = format!("{} {} {}", f[3], f[4], f[9..].join(" "));
        rows.insert(f[2].to_string(), payload);
    }
    (samples, rows)
}
