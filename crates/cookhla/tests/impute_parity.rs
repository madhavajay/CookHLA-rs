//! End-to-end parity for the IMPUTE→CONSENSUS back half, all in Rust: drive `beagle-rs` for the
//! nine imputations via [`cookhla::impute`] and check the resulting HLA calls against the oracle
//! golden. Skips loudly if the golden CONVERT_IN inputs or the `beagle-rs` binary are absent.

use std::path::PathBuf;

use cookhla::impute::{default_beagle_bin, impute_and_call, ExonInputs, ImputeConfig};

fn golden_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/example/golden")
        .canonicalize()
        .unwrap_or_else(|_| PathBuf::from("/nonexistent"))
}

#[test]
fn beagle_rs_back_half_reproduces_golden_calls() {
    let g = golden_dir();
    let gt = g.join("1958BC+HM_CEU_REF.MHC.QC.vcf");
    let alleles = g.join("1958BC+HM_CEU_REF.MHC.HLA_IMPUTATION_OUT.alleles");
    let Some(beagle_bin) = default_beagle_bin() else {
        eprintln!("SKIP: beagle-rs binary not found. Build it: (cd repos/beagle-rs && cargo build --release -p beagle-rs-cli)");
        return;
    };
    if !gt.exists() || !alleles.exists() {
        eprintln!(
            "SKIP: golden CONVERT_IN inputs not found at {}. Run docker/oracle-run.sh.",
            g.display()
        );
        return;
    }

    let err = std::fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../repos/CookHLA/example/AGM.1958BC+HM_CEU_REF.aver.erate"),
    )
    .unwrap()
    .split_whitespace()
    .next()
    .unwrap()
    .to_string();

    let exons = (2u8..=4)
        .map(|e| ExonInputs {
            exon: e,
            ref_phased_vcf: g.join(format!("HM_CEU_REF.exon{e}.phased.vcf")),
            map: g.join(format!(
                "AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB.exon{e}.txt"
            )),
        })
        .collect();

    let cfg = ImputeConfig {
        beagle_bin,
        gt_vcf: gt,
        exons,
        overlaps: vec!["0.5".into(), "1".into(), "1.5".into()],
        err,
        window: 5.0,
        ne: 10000,
        nthreads: 1,
        seed: 99999,
        workdir: std::env::temp_dir().join("cookhla_impute_parity"),
    };

    let calls = impute_and_call(&cfg).expect("impute_and_call");

    // Golden calls: cols FID _ gene 2d 4d.
    let golden = std::fs::read_to_string(&alleles).unwrap();
    let golden_calls: Vec<(String, String, String, String)> = golden
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            (f[0].into(), f[2].into(), f[3].into(), f[4].into())
        })
        .collect();

    assert_eq!(calls.len(), golden_calls.len(), "row count");

    let mut mism = 0;
    for (c, g) in calls.iter().zip(&golden_calls) {
        let two = format!("{},{}", c.two_digit.0, c.two_digit.1);
        let four = format!("{},{}", c.four_digit.0, c.four_digit.1);
        if !(c.fid == g.0 && c.gene == g.1 && two == g.2 && four == g.3) {
            if mism < 10 {
                eprintln!(
                    "MISMATCH rust [{} {} {} {}] vs golden [{} {} {} {}]",
                    c.fid, c.gene, two, four, g.0, g.1, g.2, g.3
                );
            }
            mism += 1;
        }
    }
    eprintln!(
        "beagle-rs back-half: {} calls, {} mismatches",
        calls.len(),
        mism
    );
    assert_eq!(
        mism, 0,
        "beagle-rs + consensus must reproduce the golden HLA calls"
    );
}
