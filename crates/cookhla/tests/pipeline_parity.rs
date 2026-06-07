//! Full target-side pipeline parity: QC → CONVERT_IN → IMPUTE (beagle-rs ×9) → CONSENSUS, run as
//! one Rust call from the golden NoAmbig target + the example reference + the (golden) per-exon
//! panels & maps. The resulting `.alleles` HLA calls must match the oracle golden exactly.
//! Skips loudly if plink / beagle-rs / the golden fixtures are absent.

use std::path::PathBuf;

use cookhla::front::Plink;
use cookhla::impute::default_beagle_bin;
use cookhla::pipeline::{
    run, run_building_panels, run_from_raw, ExonPanel, FullPipelineInputs, PipelineInputs,
};

fn read_err() -> String {
    std::fs::read_to_string(
        manifest().join("../../repos/CookHLA/example/AGM.1958BC+HM_CEU_REF.aver.erate"),
    )
    .unwrap()
    .split_whitespace()
    .next()
    .unwrap()
    .to_string()
}

fn golden_calls() -> Vec<(String, String, String, String)> {
    let p = golden_dir().join("1958BC+HM_CEU_REF.MHC.HLA_IMPUTATION_OUT.alleles");
    std::fs::read_to_string(p)
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            (f[0].into(), f[2].into(), f[3].into(), f[4].into())
        })
        .collect()
}

fn count_mismatches(calls: &[cookhla::consensus::AlleleCall]) -> usize {
    let golden = golden_calls();
    assert_eq!(calls.len(), golden.len(), "call count");
    let mut mism = 0;
    for (c, gc) in calls.iter().zip(&golden) {
        let two = format!("{},{}", c.two_digit.0, c.two_digit.1);
        let four = format!("{},{}", c.four_digit.0, c.four_digit.1);
        if !(c.fid == gc.0 && c.gene == gc.1 && two == gc.2 && four == gc.3) {
            if mism < 10 {
                eprintln!(
                    "MISMATCH rust [{} {} {} {}] vs golden [{} {} {} {}]",
                    c.fid, c.gene, two, four, gc.0, gc.1, gc.2, gc.3
                );
            }
            mism += 1;
        }
    }
    mism
}

fn manifest() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}
fn golden_dir() -> PathBuf {
    manifest()
        .join("../../fixtures/example/golden")
        .canonicalize()
        .unwrap_or_else(|_| PathBuf::from("/nonexistent"))
}

#[test]
fn full_target_pipeline_reproduces_golden_calls() {
    let g = golden_dir();
    let target = g.join("1958BC.hg19.COPY.LiftDown_hg18.NoAmbig");
    let golden_alleles = g.join("1958BC+HM_CEU_REF.MHC.HLA_IMPUTATION_OUT.alleles");
    let (Some(plink), Some(beagle_bin)) = (Plink::locate(), default_beagle_bin()) else {
        eprintln!("SKIP: plink and/or beagle-rs not found.");
        return;
    };
    if !target.with_extension("bed").exists() || !golden_alleles.exists() {
        eprintln!(
            "SKIP: golden fixtures not at {}. Run docker/oracle-run.sh.",
            g.display()
        );
        return;
    }

    let example_ref = manifest().join("../../repos/CookHLA/example/HM_CEU_REF");
    let err = std::fs::read_to_string(
        manifest().join("../../repos/CookHLA/example/AGM.1958BC+HM_CEU_REF.aver.erate"),
    )
    .unwrap()
    .split_whitespace()
    .next()
    .unwrap()
    .to_string();

    let panels = (2u8..=4)
        .map(|e| ExonPanel {
            exon: e,
            phased_vcf: g.join(format!("HM_CEU_REF.exon{e}.phased.vcf")),
            map: g.join(format!(
                "AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB.exon{e}.txt"
            )),
        })
        .collect();

    let workdir = std::env::temp_dir().join("cookhla_pipeline_parity");
    let inputs = PipelineInputs {
        plink,
        beagle_bin,
        target_prefix: target,
        reference_prefix: example_ref.clone(),
        reference_markers: example_ref.with_extension("markers"),
        panels,
        overlaps: vec!["0.5".into(), "1".into(), "1.5".into()],
        err,
        window: 5.0,
        ne: 10000,
        nthreads: 1,
        seed: 99999,
        small_sample: true,
        workdir: workdir.clone(),
        out_prefix: workdir.join("out"),
    };

    let calls = run(&inputs).expect("pipeline run");

    // Compare HLA calls to the golden .alleles.
    let golden = std::fs::read_to_string(&golden_alleles).unwrap();
    let golden_calls: Vec<(String, String, String, String)> = golden
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            (f[0].into(), f[2].into(), f[3].into(), f[4].into())
        })
        .collect();

    assert_eq!(calls.len(), golden_calls.len(), "call count");
    let mut mism = 0;
    for (c, gc) in calls.iter().zip(&golden_calls) {
        let two = format!("{},{}", c.two_digit.0, c.two_digit.1);
        let four = format!("{},{}", c.four_digit.0, c.four_digit.1);
        if !(c.fid == gc.0 && c.gene == gc.1 && two == gc.2 && four == gc.3) {
            if mism < 10 {
                eprintln!(
                    "MISMATCH rust [{} {} {} {}] vs golden [{} {} {} {}]",
                    c.fid, c.gene, two, four, gc.0, gc.1, gc.2, gc.3
                );
            }
            mism += 1;
        }
    }
    eprintln!(
        "FULL PIPELINE (QC→CONVERT_IN→IMPUTE→CONSENSUS): {} calls, {} mismatches",
        calls.len(),
        mism
    );
    assert_eq!(
        mism, 0,
        "full target pipeline must reproduce the golden HLA calls"
    );
}

#[test]
fn end_to_end_building_panels_and_maps_reproduces_golden_calls() {
    // Like the above, but ALSO builds the exon panels + per-exon genetic maps natively from the
    // reference + the provided AGM (no precomputed panels). The only inputs are the NoAmbig
    // target, the reference panel, and the provided whole-region AGM — everything else is built.
    let g = golden_dir();
    let target = g.join("1958BC.hg19.COPY.LiftDown_hg18.NoAmbig");
    let (Some(plink), Some(beagle_bin)) = (Plink::locate(), default_beagle_bin()) else {
        eprintln!("SKIP: plink and/or beagle-rs not found.");
        return;
    };
    if !target.with_extension("bed").exists()
        || !g
            .join("1958BC+HM_CEU_REF.MHC.HLA_IMPUTATION_OUT.alleles")
            .exists()
    {
        eprintln!("SKIP: golden fixtures absent. Run docker/oracle-run.sh.");
        return;
    }

    let example_ref = manifest().join("../../repos/CookHLA/example/HM_CEU_REF");
    let provided_agm =
        manifest().join("../../repos/CookHLA/example/AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB");
    let workdir = std::env::temp_dir().join("cookhla_e2e_full");

    let inputs = FullPipelineInputs {
        plink,
        beagle_bin,
        target_prefix: target,
        reference_prefix: example_ref,
        provided_agm,
        overlaps: vec!["0.5".into(), "1".into(), "1.5".into()],
        err: read_err(),
        window: 5.0,
        ne: 10000,
        nthreads: 1,
        seed: 99999,
        small_sample: true,
        workdir: workdir.clone(),
        out_prefix: workdir.join("out"),
    };

    let calls = run_building_panels(&inputs).expect("run_building_panels");
    let mism = count_mismatches(&calls);
    eprintln!(
        "END-TO-END (build panels+maps → QC→CONVERT_IN→IMPUTE→CONSENSUS): {} calls, {} mismatches",
        calls.len(),
        mism
    );
    assert_eq!(
        mism, 0,
        "native end-to-end must reproduce the golden HLA calls"
    );
}

#[test]
fn raw_end_to_end_from_hg19_input_reproduces_golden_calls() {
    // The COMPLETE default pipeline from the RAW hg19 target: FixInput (liftover + subset +
    // relabel + de-ambiguate) → build panels & maps → QC → CONVERT_IN → IMPUTE → CONSENSUS.
    // Inputs: the raw target PLINK files, the reference panel, and the provided AGM. Nothing
    // precomputed. Must reproduce the oracle's golden HLA calls.
    let g = golden_dir();
    let (Some(plink), Some(beagle_bin)) = (Plink::locate(), default_beagle_bin()) else {
        eprintln!("SKIP: plink and/or beagle-rs not found.");
        return;
    };
    if !g
        .join("1958BC+HM_CEU_REF.MHC.HLA_IMPUTATION_OUT.alleles")
        .exists()
    {
        eprintln!("SKIP: golden fixtures absent. Run docker/oracle-run.sh.");
        return;
    }

    let raw_target = manifest().join("../../repos/CookHLA/example/1958BC.hg19");
    let example_ref = manifest().join("../../repos/CookHLA/example/HM_CEU_REF");
    let provided_agm =
        manifest().join("../../repos/CookHLA/example/AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB");
    let workdir = std::env::temp_dir().join("cookhla_raw_e2e");

    let inputs = FullPipelineInputs {
        plink,
        beagle_bin,
        target_prefix: raw_target,
        reference_prefix: example_ref,
        provided_agm,
        overlaps: vec!["0.5".into(), "1".into(), "1.5".into()],
        err: read_err(),
        window: 5.0,
        ne: 10000,
        nthreads: 1,
        seed: 99999,
        small_sample: true,
        workdir: workdir.clone(),
        out_prefix: workdir.join("out"),
    };

    let calls = run_from_raw(&inputs, "19").expect("run_from_raw");
    let mism = count_mismatches(&calls);
    eprintln!(
        "RAW END-TO-END (FixInput → panels → QC → CONVERT_IN → IMPUTE → CONSENSUS): {} calls, {} mismatches",
        calls.len(),
        mism
    );
    assert_eq!(
        mism, 0,
        "raw end-to-end must reproduce the golden HLA calls"
    );
}
