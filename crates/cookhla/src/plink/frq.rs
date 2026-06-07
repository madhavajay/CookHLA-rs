//! PLINK `.frq` — allele frequencies from `plink --freq`. A header row then space-padded
//! columns: `CHR  SNP  A1  A2  MAF  NCHROBS`.
//!
//! CookHLA renames the header tokens with `sed` (`A1`→`A1I`, `MAF`→`MAF_I`) and then joins the
//! target and reference `.frq` on the `SNP` column via `merge_tables.pl`, comparing alleles and
//! computing `1 - MAF` for flips. We keep `MAF` as text to preserve PLINK's exact rendering
//! (e.g. `0.04545`) for byte-comparable output, and expose a typed accessor for the math.

use std::path::Path;

use anyhow::{Context, Result};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FrqRecord {
    pub chr: String,
    pub snp: String,
    pub a1: String,
    pub a2: String,
    pub maf: String,
    pub nchrobs: String,
}

impl FrqRecord {
    pub fn parse(line: &str) -> Result<Option<Self>> {
        let mut f = crate::ws_fields(line);
        let Some(chr) = f.next() else {
            return Ok(None);
        };
        let snp = f.next().context("frq: missing SNP (col 2)")?;
        let a1 = f.next().context("frq: missing A1 (col 3)")?;
        let a2 = f.next().context("frq: missing A2 (col 4)")?;
        let maf = f.next().context("frq: missing MAF (col 5)")?;
        let nchrobs = f.next().context("frq: missing NCHROBS (col 6)")?;
        Ok(Some(FrqRecord {
            chr: chr.to_owned(),
            snp: snp.to_owned(),
            a1: a1.to_owned(),
            a2: a2.to_owned(),
            maf: maf.to_owned(),
            nchrobs: nchrobs.to_owned(),
        }))
    }

    /// Minor-allele frequency as a number (for flip/strand math).
    pub fn maf_f64(&self) -> Result<f64> {
        self.maf
            .parse()
            .with_context(|| format!("frq: bad MAF {:?} for {}", self.maf, self.snp))
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Frq {
    pub records: Vec<FrqRecord>,
}

impl Frq {
    /// Parse a `.frq`, skipping the `CHR SNP ...` header line if present.
    pub fn parse(text: &str) -> Result<Self> {
        let mut records = Vec::new();
        for (i, line) in text.lines().enumerate() {
            let trimmed = line.trim_start();
            if trimmed.starts_with("CHR") {
                continue; // header
            }
            if let Some(r) = FrqRecord::parse(line)
                .with_context(|| format!("frq: parse error on line {}", i + 1))?
            {
                records.push(r);
            }
        }
        Ok(Frq { records })
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("frq: reading {}", path.display()))?;
        Self::parse(&text).with_context(|| format!("frq: parsing {}", path.display()))
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    // Real header + rows from repos/CookHLA/example/HM_CEU_REF.FRQ.frq (space-padded).
    const SAMPLE: &str = " CHR                       SNP   A1   A2          MAF  NCHROBS\n\
        \x20  6                 rs7772982    G    A        0.186      172\n\
        \x20  6                 rs9380122    C    T      0.04545      176\n";

    #[test]
    fn skips_header_and_parses_padded_columns() {
        let frq = Frq::parse(SAMPLE).unwrap();
        assert_eq!(frq.records.len(), 2);
        let r = &frq.records[0];
        assert_eq!(r.snp, "rs7772982");
        assert_eq!(r.a1, "G");
        assert_eq!(r.maf, "0.186");
        assert_eq!(r.nchrobs, "172");
    }

    #[test]
    fn maf_parses_as_number() {
        let frq = Frq::parse(SAMPLE).unwrap();
        assert!((frq.records[1].maf_f64().unwrap() - 0.04545).abs() < 1e-9);
    }
}
