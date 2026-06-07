//! Pipeline orchestrator — wires the ported stages into one run:
//! QC → CONVERT_IN (target) → IMPUTE (`beagle-rs` ×9) → CONSENSUS → `.alleles`.
//!
//! The reference-side inputs (per-exon phased panels + per-exon adaptive genetic maps) are taken
//! as paths: built by [`crate::front`]'s panel stage (in progress) or, for the example, the
//! precomputed ones. This is the verified target-side chain — it reproduces the golden HLA calls
//! from QC'd genotypes + panels + maps.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::bgl::{Bgl, Markers};
use crate::consensus::AlleleCall;
use crate::front::agm::{agm_to_text, build_exon_maps};
use crate::front::fixinput::fix_input;
use crate::front::panel::{make_exon234, make_exon_panel_vcf};
use crate::front::qc::run_qc;
use crate::front::{convert_in::convert_in_target, Plink};
use crate::impute::{impute_and_call, ExonInputs, ImputeConfig};

/// One exon panel: the phased reference VCF and the refined per-exon genetic map.
#[derive(Debug, Clone)]
pub struct ExonPanel {
    pub exon: u8,
    pub phased_vcf: PathBuf,
    pub map: PathBuf,
}

/// Inputs for a target-side imputation run (panels/maps provided).
#[derive(Debug, Clone)]
pub struct PipelineInputs {
    pub plink: Plink,
    pub beagle_bin: PathBuf,
    /// De-ambiguated target prefix (`.bed/.bim/.fam`), i.e. the FixInput/NoAmbig output.
    pub target_prefix: PathBuf,
    pub reference_prefix: PathBuf,
    pub reference_markers: PathBuf,
    pub panels: Vec<ExonPanel>, // exon 2,3,4
    pub overlaps: Vec<String>,  // e.g. ["0.5","1","1.5"]
    pub err: String,
    pub window: f64,
    pub ne: u32,
    pub nthreads: u32,
    pub seed: u64,
    pub small_sample: bool,
    pub workdir: PathBuf,
    /// Output prefix; `<out>.alleles` is written.
    pub out_prefix: PathBuf,
}

/// Run QC → CONVERT_IN → IMPUTE → CONSENSUS and write `<out_prefix>.alleles`.
pub fn run(inputs: &PipelineInputs) -> Result<Vec<AlleleCall>> {
    std::fs::create_dir_all(&inputs.workdir).ok();
    let mhc_prefix = inputs.workdir.join("MHC");

    // 1. QC
    let qc = run_qc(
        &inputs.plink,
        &inputs.target_prefix,
        &inputs.reference_prefix,
        &mhc_prefix,
        inputs.small_sample,
    )
    .context("QC stage")?;

    // 2. CONVERT_IN (target) → MHC.QC.vcf
    let vcf_text = convert_in_target(&qc, &inputs.reference_markers, "6").context("CONVERT_IN")?;
    let gt_vcf = inputs.workdir.join("MHC.QC.vcf");
    std::fs::write(&gt_vcf, vcf_text).with_context(|| format!("writing {}", gt_vcf.display()))?;

    // 3. IMPUTE (beagle-rs ×9, parallel) + 4. CONSENSUS
    let exons = inputs
        .panels
        .iter()
        .map(|p| ExonInputs {
            exon: p.exon,
            ref_phased_vcf: p.phased_vcf.clone(),
            map: p.map.clone(),
        })
        .collect();
    let cfg = ImputeConfig {
        beagle_bin: inputs.beagle_bin.clone(),
        gt_vcf,
        exons,
        overlaps: inputs.overlaps.clone(),
        err: inputs.err.clone(),
        window: inputs.window,
        ne: inputs.ne,
        nthreads: inputs.nthreads,
        seed: inputs.seed,
        workdir: inputs.workdir.join("impute"),
    };
    let calls = impute_and_call(&cfg).context("IMPUTE/CONSENSUS")?;

    // 5. Write .alleles
    let alleles_path = with_suffix(&inputs.out_prefix, ".alleles");
    let mut body = String::new();
    for c in &calls {
        body.push_str(&c.to_alleles_line());
        body.push('\n');
    }
    std::fs::write(&alleles_path, body)
        .with_context(|| format!("writing {}", alleles_path.display()))?;

    Ok(calls)
}

fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

/// Inputs for a run that also **builds the reference panels and per-exon maps natively** from the
/// reference panel + the provided adaptive genetic map (no precomputed panels needed).
#[derive(Debug, Clone)]
pub struct FullPipelineInputs {
    pub plink: Plink,
    pub beagle_bin: PathBuf,
    pub target_prefix: PathBuf,
    /// Reference prefix: `.bgl.phased`, `.markers`, `.bim`, `.FRQ.frq`.
    pub reference_prefix: PathBuf,
    /// Provided adaptive genetic map (`*.mach_step.avg.clpsB`).
    pub provided_agm: PathBuf,
    pub overlaps: Vec<String>,
    pub err: String,
    pub window: f64,
    pub ne: u32,
    pub nthreads: u32,
    pub seed: u64,
    pub small_sample: bool,
    pub workdir: PathBuf,
    pub out_prefix: PathBuf,
}

/// The complete default pipeline from a **raw** target: FixInput (liftover hg→18 + subset +
/// relabel + de-ambiguate) → build panels & maps → QC → CONVERT_IN → IMPUTE → CONSENSUS. The
/// provided adaptive genetic map is used (MakeGeneticMap/MACH is not run). `target_hg` is the
/// target build (`"18"`, `"19"`, `"38"`).
pub fn run_from_raw(inputs: &FullPipelineInputs, target_hg: &str) -> Result<Vec<AlleleCall>> {
    std::fs::create_dir_all(&inputs.workdir).ok();

    // FixInput: raw target → de-ambiguated, reference-matched target.
    let noambig = fix_input(
        &inputs.plink,
        &inputs.target_prefix,
        target_hg,
        &inputs.reference_prefix,
        &inputs.workdir.join("target.COPY"),
    )
    .context("FixInput")?;

    let full = FullPipelineInputs {
        target_prefix: noambig,
        ..inputs.clone()
    };
    run_building_panels(&full)
}

/// Build the exon panels + per-exon maps from the reference, then run the target-side pipeline.
/// This is everything except FixInput (liftover) and MakeGeneticMap (the provided AGM is used).
pub fn run_building_panels(inputs: &FullPipelineInputs) -> Result<Vec<AlleleCall>> {
    std::fs::create_dir_all(&inputs.workdir).ok();

    // --- build exon panels (reference side) ---
    let ref_bgl = Bgl::read(with_suffix(&inputs.reference_prefix, ".bgl.phased"))
        .context("reading reference .bgl.phased")?;
    let ref_markers = Markers::read(with_suffix(&inputs.reference_prefix, ".markers"))
        .context("reading reference .markers")?;
    let (e234_bgl, e234_markers) = make_exon234(&ref_bgl, &ref_markers);

    // --- build per-exon adaptive genetic maps ---
    let maps =
        build_exon_maps(&inputs.provided_agm, &e234_markers).context("building exon maps")?;

    let mut panels = Vec::new();
    for exon in [2u8, 3, 4] {
        let vcf = make_exon_panel_vcf(exon, &e234_bgl, &e234_markers)?;
        let vcf_path = inputs.workdir.join(format!("ref.exon{exon}.phased.vcf"));
        std::fs::write(&vcf_path, vcf)
            .with_context(|| format!("writing {}", vcf_path.display()))?;
        let map_path = inputs.workdir.join(format!("ref.exon{exon}.map.txt"));
        std::fs::write(&map_path, agm_to_text(&maps[&exon]))
            .with_context(|| format!("writing {}", map_path.display()))?;
        panels.push(ExonPanel {
            exon,
            phased_vcf: vcf_path,
            map: map_path,
        });
    }

    run(&PipelineInputs {
        plink: inputs.plink.clone(),
        beagle_bin: inputs.beagle_bin.clone(),
        target_prefix: inputs.target_prefix.clone(),
        reference_prefix: inputs.reference_prefix.clone(),
        reference_markers: with_suffix(&inputs.reference_prefix, ".markers"),
        panels,
        overlaps: inputs.overlaps.clone(),
        err: inputs.err.clone(),
        window: inputs.window,
        ne: inputs.ne,
        nthreads: inputs.nthreads,
        seed: inputs.seed,
        small_sample: inputs.small_sample,
        workdir: inputs.workdir.clone(),
        out_prefix: inputs.out_prefix.clone(),
    })
}
