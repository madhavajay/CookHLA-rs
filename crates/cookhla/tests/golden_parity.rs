//! End-to-end parity for the consensus caller against the **real** CookHLA golden output.
//!
//! Reproduces `CONVERT_OUT` from the nine real Beagle imputation VCFs captured by
//! `docker/oracle-run.sh` (→ `fixtures/example/golden/`) and checks the result against the
//! reference `.alleles`: HLA *calls* must match exactly; posterior probabilities within epsilon.
//!
//! The fixtures are large and git-ignored, so this test **skips loudly** (samtools-rs honest-gate
//! lesson) when they are absent rather than silently passing. Regenerate them with
//! `docker build -f docker/Dockerfile.oracle -t cookhla-oracle . && docker/oracle-run.sh`.

use std::collections::BTreeMap;
use std::path::PathBuf;

use cookhla::consensus::{convert_out, order_nine, ImputationVcf};

fn golden_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/example/golden")
        .canonicalize()
        .unwrap_or_else(|_| PathBuf::from("/nonexistent"))
}

/// Parse `(exon, overlap)` out of a name like `...exon2.0.5.raw_imputation_out.vcf` — note the
/// overlap can itself contain a dot (`0.5`, `1.5`), so split on the literal markers, not `.`.
fn parse_exon_overlap(name: &str) -> Option<(u8, String)> {
    let idx = name.find("exon")?;
    let rest = &name[idx + 4..]; // e.g. "2.0.5.raw_imputation_out.vcf"
    let dot = rest.find('.')?;
    let exon: u8 = rest[..dot].parse().ok()?;
    let rest2 = &rest[dot + 1..]; // "0.5.raw_imputation_out.vcf"
    let ov_end = rest2.find(".raw_imputation_out")?;
    Some((exon, rest2[..ov_end].to_string()))
}

/// Columns of a `.alleles` line: FID IID gene 2d 4d pp1 pp2 conf.
struct AllelesLine {
    fid: String,
    gene: String,
    two_digit: String,
    four_digit: String,
    pp: [f64; 3],
}

fn parse_alleles(text: &str) -> Vec<AllelesLine> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            AllelesLine {
                fid: f[0].to_string(),
                gene: f[2].to_string(),
                two_digit: f[3].to_string(),
                four_digit: f[4].to_string(),
                pp: [
                    f[5].parse().unwrap(),
                    f[6].parse().unwrap(),
                    f[7].parse().unwrap(),
                ],
            }
        })
        .collect()
}

#[test]
fn consensus_matches_golden_alleles() {
    let dir = golden_dir();
    let alleles_path = dir.join("1958BC+HM_CEU_REF.MHC.HLA_IMPUTATION_OUT.alleles");
    if !alleles_path.exists() {
        eprintln!(
            "SKIP consensus_matches_golden_alleles: golden not found at {}.\n\
             Regenerate: docker build -f docker/Dockerfile.oracle -t cookhla-oracle . && docker/oracle-run.sh",
            alleles_path.display()
        );
        return;
    }

    // Imputation VCFs come from COOKHLA_IMP_DIR if set (e.g. a beagle-rs run), else the golden
    // dir (the oracle's own Beagle output). The golden `.alleles` is always the oracle's.
    let imp_dir = std::env::var("COOKHLA_IMP_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| dir.clone());

    // Load the nine imputation VCFs, keyed by (exon, overlap).
    let mut by_key: BTreeMap<(u8, String), ImputationVcf> = BTreeMap::new();
    for entry in std::fs::read_dir(&imp_dir).unwrap() {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        if !name.ends_with(".raw_imputation_out.vcf") {
            continue;
        }
        let Some(key) = parse_exon_overlap(&name) else {
            continue;
        };
        let text = std::fs::read_to_string(&path).unwrap();
        by_key.insert(key, ImputationVcf::parse(&text).unwrap());
    }
    assert_eq!(
        by_key.len(),
        9,
        "expected 9 imputation VCFs, found {}",
        by_key.len()
    );

    let nine = order_nine(by_key).unwrap();
    let calls = convert_out(&nine);

    let golden = parse_alleles(&std::fs::read_to_string(&alleles_path).unwrap());

    assert_eq!(
        calls.len(),
        golden.len(),
        "row count: rust {} vs golden {}",
        calls.len(),
        golden.len()
    );

    let mut call_mismatches = 0;
    let mut max_prob_drift = 0.0f64;
    for (c, g) in calls.iter().zip(&golden) {
        let c_two = format!("{},{}", c.two_digit.0, c.two_digit.1);
        let c_four = format!("{},{}", c.four_digit.0, c.four_digit.1);
        let calls_match =
            c.fid == g.fid && c.gene == g.gene && c_two == g.two_digit && c_four == g.four_digit;
        if !calls_match {
            if call_mismatches < 15 {
                eprintln!(
                    "CALL MISMATCH: rust [{} {} {} {}] vs golden [{} {} {} {}]",
                    c.fid, c.gene, c_two, c_four, g.fid, g.gene, g.two_digit, g.four_digit
                );
            }
            call_mismatches += 1;
        }
        for (cp, gp) in [c.pp1, c.pp2, c.conf].iter().zip(&g.pp) {
            max_prob_drift = max_prob_drift.max((cp - gp).abs());
        }
    }

    // Two modes:
    //  - default (golden's own Beagle 5.1 VCFs): pure consensus parity — probabilities must match
    //    the R original to floating-point precision.
    //  - COOKHLA_IMP_DIR set (e.g. beagle-rs 5.5 output): end-to-end engine parity — HLA *calls*
    //    must match exactly; probabilities legitimately drift with the Beagle version.
    let external_engine = std::env::var("COOKHLA_IMP_DIR").is_ok();
    eprintln!(
        "consensus parity ({}): {} calls, {} call-mismatches, max |Δprob| = {:.3e}",
        if external_engine {
            "engine: beagle-rs vs golden"
        } else {
            "consensus vs R"
        },
        calls.len(),
        call_mismatches,
        max_prob_drift
    );
    assert_eq!(
        call_mismatches, 0,
        "HLA calls must match the golden exactly"
    );
    if external_engine {
        // Beagle 5.5 vs the oracle's 5.1: a few % posterior drift is expected; calls must hold.
        assert!(
            max_prob_drift < 0.2,
            "probability drift {max_prob_drift:.3e} unexpectedly large for a version gap"
        );
    } else {
        assert!(
            max_prob_drift < 1e-3,
            "consensus probability drift {max_prob_drift:.3e} exceeds 1e-3 vs the R original"
        );
    }
}
