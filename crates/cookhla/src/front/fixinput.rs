//! FixInput (`checkInput.py::FixInput` + `exclude_Ambiguous_SNP`) — prepare the raw target for
//! imputation: lift it down to hg18, subset to the reference markers by base position, relabel to
//! the reference marker names, and drop ambiguous (`{A,T}`/`{G,C}`) SNPs.
//!
//! `plink` does the genotype-matrix subset/relabel; the liftover ([`super::liftover`]) and the
//! base-position matching are native Rust. The target alleles are kept (only the marker *name* and
//! *position* are updated to the reference); the rare 0-allele fix-up in `UpdateInput` is applied
//! when present.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use super::liftover::LiftOver;
use super::Plink;
use crate::plink::Bim;

const ACGT: [&str; 4] = ["A", "C", "G", "T"];

fn is_mkref(label: &str) -> bool {
    label.starts_with("AA_")
        || label.starts_with("HLA_")
        || label.starts_with("SNP_")
        || label.starts_with("INS_")
}

fn is_ambiguous(a: &str, b: &str) -> bool {
    matches!((a, b), ("A", "T") | ("T", "A") | ("C", "G") | ("G", "C"))
}

/// Run FixInput. `hg` is the target build (`"18"`, `"19"`, `"38"`); for `"19"`/`"38"` the target is
/// lifted to hg18 first (hg38 → hg19 → hg18 — the second hop is not yet vendored). Returns the
/// `NoAmbig` prefix (`.bed/.bim/.fam`).
pub fn fix_input(
    plink: &Plink,
    target_prefix: &Path,
    hg: &str,
    reference_prefix: &Path,
    out_prefix: &Path,
) -> Result<PathBuf> {
    let workdir = out_prefix.parent().unwrap_or(Path::new(".")).to_path_buf();
    std::fs::create_dir_all(&workdir).ok();

    // --- 1. Lift the target .bim down to hg18 (or pass through if already hg18) ---
    let target_bim = Bim::read(with(target_prefix, ".bim"))?;
    let lo = (hg != "18").then(LiftOver::new);
    let mut lifted = target_bim.clone();
    for r in &mut lifted.records {
        if let Some(lo) = &lo {
            match lo.convert(r.bp) {
                Some(p) => r.bp = p,
                None => r.bp = -1, // failed lift → sentinel; excluded below
            }
        }
    }
    let lifted_bim_path = with(out_prefix, ".lifted.bim");
    lifted.write(&lifted_bim_path)?;

    // --- 2. Match target → reference by hg18 base position (reference SNPs only) ---
    let ref_bim = Bim::read(with(reference_prefix, ".bim"))?;
    let mut ref_by_bp: HashMap<i64, &str> = HashMap::new();
    for r in &ref_bim.records {
        if is_mkref(&r.id) {
            continue; // exclude MakeReference markers (AA_/HLA_/SNP_/INS_)
        }
        ref_by_bp.entry(r.bp).or_insert(r.id.as_str()); // first per BP (drop_duplicates)
    }

    // Build the extract list (target labels) + update_name (target → reference), one per BP.
    let mut extract = String::new();
    let mut update_name = String::new();
    let mut seen_bp: HashSet<i64> = HashSet::new();
    for r in &lifted.records {
        if r.bp < 0 {
            continue;
        }
        if let Some(&ref_label) = ref_by_bp.get(&r.bp) {
            if !seen_bp.insert(r.bp) {
                continue; // one target marker per BP
            }
            extract.push_str(&r.id);
            extract.push('\n');
            update_name.push_str(&format!("{}\t{}\n", r.id, ref_label));
        }
    }
    let extract_path = with(out_prefix, ".extract");
    let update_name_path = with(out_prefix, ".update_name");
    std::fs::write(&extract_path, extract)?;
    std::fs::write(&update_name_path, update_name)?;

    // --- 3. plink: subset to matched markers, then relabel to reference names ---
    let subset = with(out_prefix, ".subset");
    let subset_s = subset.to_string_lossy().to_string();
    plink
        .run(&[
            "--bed",
            &with(target_prefix, ".bed").to_string_lossy(),
            "--bim",
            &lifted_bim_path.to_string_lossy(),
            "--fam",
            &with(target_prefix, ".fam").to_string_lossy(),
            "--extract",
            &extract_path.to_string_lossy(),
            "--make-bed",
            "--out",
            &subset_s,
        ])
        .context("FixInput: subset by base position")?;

    let fixed = out_prefix.to_string_lossy().to_string();
    plink
        .run(&[
            "--bfile",
            &subset_s,
            "--update-name",
            &update_name_path.to_string_lossy(),
            "--make-bed",
            "--out",
            &fixed,
        ])
        .context("FixInput: relabel to reference marker names")?;

    // --- 4. Exclude ambiguous SNPs ({A,T}/{G,C}) → NoAmbig ---
    let fixed_bim = Bim::read(with(out_prefix, ".bim"))?;
    let mut ambig = String::new();
    for r in &fixed_bim.records {
        let (a1, a2) = (r.a1.as_str(), r.a2.as_str());
        if ACGT.contains(&a1) && ACGT.contains(&a2) && is_ambiguous(a1, a2) {
            ambig.push_str(&r.id);
            ambig.push('\n');
        }
    }
    let ambig_path = with(out_prefix, ".ambig");
    std::fs::write(&ambig_path, ambig)?;

    let noambig = with(out_prefix, ".NoAmbig");
    let noambig_s = noambig.to_string_lossy().to_string();
    plink
        .run(&[
            "--bfile",
            &fixed,
            "--exclude",
            &ambig_path.to_string_lossy(),
            "--make-bed",
            "--keep-allele-order",
            "--out",
            &noambig_s,
        ])
        .context("FixInput: exclude ambiguous SNPs")?;

    Ok(noambig)
}

/// `prefix` + `suffix` as a path (suffix includes the dot).
fn with(prefix: &Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mkref_and_ambiguity_classification() {
        assert!(is_mkref("HLA_A_0101"));
        assert!(is_mkref("AA_A_-15"));
        assert!(!is_mkref("rs123"));
        assert!(is_ambiguous("A", "T"));
        assert!(!is_ambiguous("A", "G"));
    }
}
