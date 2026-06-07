//! Thin wrapper over the `plink` (v1.9) executable — the fast compiled tool CookHLA relies on
//! for genotype-matrix operations (`--make-bed`, `--recode`, `--freq`, `--flip`, `--extract`,
//! `--exclude`, `--geno`, region filters). We keep it as an external binary (per the staged
//! strategy) and reproduce only the text glue around it in Rust.

use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{bail, Context, Result};

/// A located `plink` binary, invoked with CookHLA's standard flags.
#[derive(Debug, Clone)]
pub struct Plink {
    bin: PathBuf,
}

impl Plink {
    pub fn new(bin: impl Into<PathBuf>) -> Self {
        Plink { bin: bin.into() }
    }

    /// Locate `plink`: `$PLINK`, else `PATH`, else the vendored `repos/CookHLA/dependency/plink`.
    pub fn locate() -> Option<Self> {
        if let Ok(p) = std::env::var("PLINK") {
            let p = PathBuf::from(p);
            if p.exists() {
                return Some(Plink::new(p));
            }
        }
        // Vendored copy relative to the crate.
        let vendored =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../repos/CookHLA/dependency/plink");
        if let Ok(p) = vendored.canonicalize() {
            return Some(Plink::new(p));
        }
        // PATH lookup.
        if let Ok(out) = Command::new("which").arg("plink").output() {
            if out.status.success() {
                let p = String::from_utf8_lossy(&out.stdout).trim().to_string();
                if !p.is_empty() {
                    return Some(Plink::new(p));
                }
            }
        }
        None
    }

    pub fn path(&self) -> &Path {
        &self.bin
    }

    /// Run `plink --noweb --silent --allow-no-sex <args...>`, failing loudly on a nonzero exit.
    pub fn run(&self, args: &[&str]) -> Result<()> {
        let status = Command::new(&self.bin)
            .args(["--noweb", "--silent", "--allow-no-sex"])
            .args(args)
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .with_context(|| format!("failed to launch plink ({})", self.bin.display()))?;
        if !status.success() {
            bail!(
                "plink failed (exit {:?}) for args {:?}",
                status.code(),
                args
            );
        }
        Ok(())
    }
}
