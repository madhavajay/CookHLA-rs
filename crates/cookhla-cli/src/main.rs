//! `cookhla` — the CookHLA-rs command-line driver.
//!
//! Argv-compatible with the reference `CookHLA.py` for the default path. The reference uses
//! argparse-style single-dash multi-character flags (`-hg`, `-ref`, `-gm`, ...); `clap` models
//! those as long options, so we rewrite the legacy spellings to `--long` form before parsing.
//!
//! Runs the complete default pipeline ([`cookhla::pipeline::run_from_raw`]): FixInput → build
//! panels & maps → QC → CONVERT_IN → IMPUTE (`beagle-rs`) → CONSENSUS → `<out>.alleles`. Requires
//! a precomputed adaptive genetic map (`-gm`/`-ae`); auto-generating it (MakeGeneticMap/MACH) and
//! the legacy Beagle 4.1 path are later milestones.

use anyhow::{Context, Result};
use clap::Parser;

/// Fast, accurate HLA imputation — a Rust port of CookHLA.
#[derive(Debug, Parser)]
#[command(name = "cookhla", version, about)]
struct Args {
    /// Common prefix of the target PLINK binary input (`-i`).
    #[arg(long, short = 'i')]
    input: String,

    /// Human-genome build of the TARGET data: 18, 19, or 38 (`-hg`).
    #[arg(long = "human-genome", value_parser = ["18", "19", "38"])]
    human_genome: String,

    /// Prefix of the SNP2HLA-formatted reference panel (`-ref`).
    #[arg(long)]
    reference: String,

    /// Output prefix (`-o`).
    #[arg(long, short = 'o')]
    out: String,

    /// Adaptive genetic-map file (`-gm`). Generated via MakeGeneticMap if omitted.
    #[arg(long = "genetic-map")]
    genetic_map: Option<String>,

    /// Average error-rate file (`-ae`).
    #[arg(long = "average-erate")]
    average_erate: Option<String>,

    /// Answer file for accuracy measurement (`-an`).
    #[arg(long)]
    answer: Option<String>,

    /// Cores for multiprocessing the 9 exon×overlap imputations (`-mp`).
    #[arg(long, default_value_t = 1)]
    multiprocess: u8,

    /// Heap memory per Beagle run, e.g. `2g` (`-mem`). Retained for CLI parity.
    #[arg(long = "java-memory", default_value = "2g")]
    java_memory: String,

    /// Threads per Beagle run (`-nth`).
    #[arg(long, default_value_t = 1)]
    nthreads: u32,

    /// Use Beagle 4.1 instead of 5.x (`-bgl4`). Legacy path (later milestone).
    #[arg(long, default_value_t = false)]
    beagle4: bool,

    /// Three overlap values in cM for Beagle 5.x (`-ol`).
    #[arg(long, num_args = 3, value_delimiter = ' ', default_values_t = [0.5, 1.0, 1.5])]
    overlap: Vec<f64>,

    /// Window size in cM for Beagle 5.x (`-w`).
    #[arg(long, short = 'w', default_value_t = 5.0)]
    window: f64,

    /// Effective population size for Beagle 5.x (`-ne`).
    #[arg(long = "effective-population-size", default_value_t = 10000)]
    ne: u32,
}

/// Rewrite the reference's single-dash multi-char flags into the `--long` spellings clap knows.
fn normalize_legacy_flags(argv: impl Iterator<Item = String>) -> Vec<String> {
    argv.map(|a| match a.as_str() {
        "-hg" => "--human-genome".to_owned(),
        "-ref" => "--reference".to_owned(),
        "-gm" => "--genetic-map".to_owned(),
        "-ae" => "--average-erate".to_owned(),
        "-an" => "--answer".to_owned(),
        "-mp" => "--multiprocess".to_owned(),
        "-mem" => "--java-memory".to_owned(),
        "-nth" => "--nthreads".to_owned(),
        "-bgl4" => "--beagle4".to_owned(),
        "-ol" => "--overlap".to_owned(),
        "-ne" => "--effective-population-size".to_owned(),
        "-macc_v2" => "--measureAcc_v2".to_owned(),
        _ => a,
    })
    .collect()
}

fn main() -> Result<()> {
    let argv = normalize_legacy_flags(std::env::args());
    let args = Args::parse_from(argv);

    eprintln!("[cookhla] CookHLA-rs (v{})", env!("CARGO_PKG_VERSION"));
    eprintln!("  input        = {}", args.input);
    eprintln!("  human_genome = hg{}", args.human_genome);
    eprintln!("  reference    = {}", args.reference);
    eprintln!("  out          = {}", args.out);
    eprintln!(
        "  genetic_map  = {}",
        args.genetic_map
            .as_deref()
            .unwrap_or("<auto: MakeGeneticMap>")
    );
    eprintln!(
        "  beagle       = {}",
        if args.beagle4 {
            "4.1 (legacy)"
        } else {
            "5.x via beagle-rs"
        }
    );
    eprintln!(
        "  overlaps     = {:?}  window = {}  ne = {}  nthreads = {}",
        args.overlap, args.window, args.ne, args.nthreads
    );

    if args.beagle4 {
        anyhow::bail!("--beagle4 (legacy Beagle 4.1 path) is not yet ported; default uses Beagle 5.5 via beagle-rs.");
    }

    // The adaptive genetic map must be provided: MakeGeneticMap/MACH is a later milestone.
    let (Some(gm), Some(ae)) = (args.genetic_map.as_ref(), args.average_erate.as_ref()) else {
        anyhow::bail!(
            "this build requires a precomputed adaptive genetic map: pass -gm <...clpsB> and -ae \
             <...aver.erate>. (Auto-generating the AGM — MakeGeneticMap/MACH — is not yet ported.)"
        );
    };

    let plink = cookhla::front::Plink::locate()
        .context("plink not found — set $PLINK or vendor repos/CookHLA/dependency/plink")?;
    let beagle_bin = cookhla::impute::default_beagle_bin()
        .context("beagle-rs not found — set $BEAGLE_RS or build repos/beagle-rs (cargo build --release -p beagle-rs-cli)")?;

    let err = std::fs::read_to_string(ae)
        .with_context(|| format!("reading --average-erate {ae}"))?
        .split_whitespace()
        .next()
        .context("--average-erate file is empty")?
        .to_string();

    // Small-sample mode when the target has < 100 samples.
    let n_samples = cookhla::plink::Fam::read(format!("{}.fam", args.input))
        .map(|f| f.n_samples())
        .unwrap_or(0);
    let small_sample = n_samples < 100;
    eprintln!(
        "  samples      = {n_samples}{}",
        if small_sample {
            " (small-sample mode)"
        } else {
            ""
        }
    );

    let out_path = std::path::PathBuf::from(&args.out);
    let workdir = out_path
        .parent()
        .map(|p| p.join(format!(".{}_work", file_stem(&args.out))))
        .unwrap_or_else(|| std::path::PathBuf::from(".cookhla_work"));

    let inputs = cookhla::pipeline::FullPipelineInputs {
        plink,
        beagle_bin,
        target_prefix: std::path::PathBuf::from(&args.input),
        reference_prefix: std::path::PathBuf::from(&args.reference),
        provided_agm: std::path::PathBuf::from(gm),
        overlaps: args.overlap.iter().map(fmt_overlap).collect(),
        err,
        window: args.window,
        ne: args.ne,
        nthreads: args.nthreads,
        seed: 99999,
        small_sample,
        workdir,
        out_prefix: out_path.clone(),
    };

    let start = std::time::Instant::now();
    let calls =
        cookhla::pipeline::run_from_raw(&inputs, &args.human_genome).context("CookHLA pipeline")?;
    let elapsed = start.elapsed();

    eprintln!(
        "[cookhla] DONE — {} HLA calls written to {}.alleles in {:.1}s",
        calls.len(),
        args.out,
        elapsed.as_secs_f64()
    );
    Ok(())
}

/// Render an overlap value the way the reference labels them (`0.5`, `1`, `1.5`).
fn fmt_overlap(x: &f64) -> String {
    if *x == x.trunc() {
        format!("{}", *x as i64)
    } else {
        format!("{x}")
    }
}

fn file_stem(prefix: &str) -> String {
    std::path::Path::new(prefix)
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "cookhla".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rewrites_legacy_flags() {
        let out = normalize_legacy_flags(
            [
                "cookhla", "-i", "t", "-hg", "19", "-ref", "r", "-o", "o", "-ne", "20000",
            ]
            .into_iter()
            .map(String::from),
        );
        assert!(out.contains(&"--reference".to_string()));
        assert!(out.contains(&"--human-genome".to_string()));
        assert!(out.contains(&"--effective-population-size".to_string()));
        // -i and -o (real short flags) are left alone.
        assert!(out.contains(&"-i".to_string()));
        assert!(out.contains(&"-o".to_string()));
    }
}
