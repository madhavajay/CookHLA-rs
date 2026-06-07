//! The IMPUTE stage (`HLA_Imputation_BEAGLE5.IMPUTE`) — run the nine `exon × overlap`
//! imputations through **`beagle-rs`** (the Rust Beagle 5.5 port) and hand the results to the
//! [`crate::consensus`] caller.
//!
//! CookHLA shells out to `java -jar beagle5.jar` nine times (JVM startup ×9, serial-ish). We
//! instead invoke the compiled `beagle-rs` binary — no JVM — and run the nine in parallel with
//! `rayon`. Verified end-to-end: from the real CONVERT_IN inputs this reproduces every golden HLA
//! call (the GP posteriors drift a few percent because `beagle-rs` is 5.5 vs the reference 5.1;
//! the *calls* are identical). In-process invocation of the `beagle-rs` library is a later
//! optimization; a native-binary subprocess is already fast and keeps the crates decoupled.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{bail, Context, Result};
use rayon::prelude::*;

use crate::consensus::{convert_out, order_nine, AlleleCall, ImputationVcf};

/// Per-exon reference inputs: the phased reference panel VCF and the refined genetic map.
#[derive(Debug, Clone)]
pub struct ExonInputs {
    pub exon: u8,
    pub ref_phased_vcf: PathBuf,
    pub map: PathBuf,
}

/// Everything needed to run the nine imputations.
#[derive(Debug, Clone)]
pub struct ImputeConfig {
    /// Path to the `beagle-rs` binary (drop-in for `beagle.jar`).
    pub beagle_bin: PathBuf,
    /// Target genotypes to impute (`gt=`), i.e. the CONVERT_IN `MHC.QC.vcf`.
    pub gt_vcf: PathBuf,
    /// The three exon panels (2, 3, 4).
    pub exons: Vec<ExonInputs>,
    /// The three overlap values (cM), e.g. `["0.5", "1", "1.5"]`.
    pub overlaps: Vec<String>,
    /// Average error rate (`err=`), the scalar from the `.aver.erate` file.
    pub err: String,
    pub window: f64,
    pub ne: u32,
    pub nthreads: u32,
    pub seed: u64,
    /// Scratch directory for the per-run output VCFs.
    pub workdir: PathBuf,
}

/// Locate a `beagle-rs` binary: `$BEAGLE_RS` if set, else the in-repo release build.
pub fn default_beagle_bin() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("BEAGLE_RS") {
        let p = PathBuf::from(p);
        if p.exists() {
            return Some(p);
        }
    }
    let candidate = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../repos/beagle-rs/target/release/beagle-rs");
    candidate.canonicalize().ok()
}

/// Run one imputation and parse its output into an [`ImputationVcf`].
fn run_one(
    cfg: &ImputeConfig,
    exon: &ExonInputs,
    overlap: &str,
) -> Result<((u8, String), ImputationVcf)> {
    let prefix = cfg
        .workdir
        .join(format!("impute.exon{}.{}", exon.exon, overlap));

    let status = Command::new(&cfg.beagle_bin)
        .arg(format!("gt={}", cfg.gt_vcf.display()))
        .arg(format!("ref={}", exon.ref_phased_vcf.display()))
        .arg(format!("out={}", prefix.display()))
        .arg("impute=true")
        .arg("gp=true")
        .arg(format!("overlap={overlap}"))
        .arg(format!("err={}", cfg.err))
        .arg(format!("map={}", exon.map.display()))
        .arg(format!("window={}", cfg.window))
        .arg(format!("ne={}", cfg.ne))
        .arg(format!("nthreads={}", cfg.nthreads))
        .arg(format!("seed={}", cfg.seed))
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .with_context(|| format!("failed to launch beagle-rs ({})", cfg.beagle_bin.display()))?;
    if !status.success() {
        bail!(
            "beagle-rs failed for exon{} overlap {} (status {:?})",
            exon.exon,
            overlap,
            status.code()
        );
    }

    // Beagle writes `<out>.vcf.gz` by literally appending — don't use `with_extension`, which
    // would mangle the dot-containing prefix (e.g. `impute.exon2.0.5`).
    let gz = PathBuf::from(format!("{}.vcf.gz", prefix.display()));
    let text = read_maybe_gzip(&gz)
        .with_context(|| format!("reading imputation output {}", gz.display()))?;
    let vcf = ImputationVcf::parse(&text)?;
    Ok(((exon.exon, overlap.to_string()), vcf))
}

/// Read a `.vcf.gz` (gzip) into text.
fn read_maybe_gzip(path: &Path) -> Result<String> {
    let bytes = std::fs::read(path)?;
    let mut out = String::new();
    flate2::read::MultiGzDecoder::new(&bytes[..])
        .read_to_string(&mut out)
        .context("gunzip imputation VCF")?;
    Ok(out)
}

/// Run all nine `exon × overlap` imputations (in parallel) and key them for the consensus.
pub fn run_nine(cfg: &ImputeConfig) -> Result<BTreeMap<(u8, String), ImputationVcf>> {
    if cfg.exons.len() != 3 || cfg.overlaps.len() != 3 {
        bail!(
            "run_nine: expected 3 exons × 3 overlaps, got {}×{}",
            cfg.exons.len(),
            cfg.overlaps.len()
        );
    }
    std::fs::create_dir_all(&cfg.workdir).ok();

    let jobs: Vec<(&ExonInputs, &String)> = cfg
        .exons
        .iter()
        .flat_map(|e| cfg.overlaps.iter().map(move |o| (e, o)))
        .collect();

    let results: Vec<Result<((u8, String), ImputationVcf)>> = jobs
        .par_iter()
        .map(|(exon, overlap)| run_one(cfg, exon, overlap))
        .collect();

    let mut map = BTreeMap::new();
    for r in results {
        let (key, vcf) = r?;
        map.insert(key, vcf);
    }
    Ok(map)
}

/// IMPUTE + CONVERT_OUT: run the nine imputations and produce the final HLA allele calls.
pub fn impute_and_call(cfg: &ImputeConfig) -> Result<Vec<AlleleCall>> {
    let nine = order_nine(run_nine(cfg)?)?;
    Ok(convert_out(&nine))
}
