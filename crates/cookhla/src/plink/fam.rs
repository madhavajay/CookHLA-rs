//! PLINK `.fam` — the sample table. Six space-separated columns:
//! `fid  iid  pid  mid  sex  pheno`.
//!
//! CookHLA mostly counts lines here (`getSampleNumbers`) and rewrites the parental columns to
//! `0` before randomizing sample order in `MakeGeneticMap` (`awk '{print $1" "$2" 0 0 "$5" "$6}'`).
//! Fields are kept as text because IDs are arbitrary tokens and `pheno` is often `-9`.

use std::fmt;
use std::path::Path;

use anyhow::{Context, Result};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FamRecord {
    pub fid: String,
    pub iid: String,
    pub pid: String,
    pub mid: String,
    pub sex: String,
    pub pheno: String,
}

impl FamRecord {
    pub fn parse(line: &str) -> Result<Option<Self>> {
        let mut f = crate::ws_fields(line);
        let Some(fid) = f.next() else {
            return Ok(None);
        };
        let iid = f.next().context("fam: missing iid (col 2)")?;
        let pid = f.next().context("fam: missing pid (col 3)")?;
        let mid = f.next().context("fam: missing mid (col 4)")?;
        let sex = f.next().context("fam: missing sex (col 5)")?;
        let pheno = f.next().context("fam: missing pheno (col 6)")?;
        Ok(Some(FamRecord {
            fid: fid.to_owned(),
            iid: iid.to_owned(),
            pid: pid.to_owned(),
            mid: mid.to_owned(),
            sex: sex.to_owned(),
            pheno: pheno.to_owned(),
        }))
    }

    /// The `MakeGeneticMap` "trick fam": zero out the parental IDs, keep the rest.
    pub fn with_parents_zeroed(&self) -> FamRecord {
        FamRecord {
            pid: "0".into(),
            mid: "0".into(),
            ..self.clone()
        }
    }
}

impl fmt::Display for FamRecord {
    /// PLINK writes `.fam` space-separated.
    fn fmt(&self, w: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            w,
            "{} {} {} {} {} {}",
            self.fid, self.iid, self.pid, self.mid, self.sex, self.pheno
        )
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Fam {
    pub records: Vec<FamRecord>,
}

impl Fam {
    pub fn parse(text: &str) -> Result<Self> {
        let mut records = Vec::new();
        for (i, line) in text.lines().enumerate() {
            if let Some(r) = FamRecord::parse(line)
                .with_context(|| format!("fam: parse error on line {}", i + 1))?
            {
                records.push(r);
            }
        }
        Ok(Fam { records })
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("fam: reading {}", path.display()))?;
        Self::parse(&text).with_context(|| format!("fam: parsing {}", path.display()))
    }

    pub fn to_text(&self) -> String {
        let mut s = String::new();
        for r in &self.records {
            s.push_str(&r.to_string());
            s.push('\n');
        }
        s
    }

    pub fn write(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        std::fs::write(path, self.to_text())
            .with_context(|| format!("fam: writing {}", path.display()))
    }

    /// `getSampleNumbers` — the line count CookHLA uses to pick small-sample mode (<100).
    pub fn n_samples(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    // Real rows from repos/CookHLA/example/HM_CEU_REF.fam.
    const SAMPLE: &str = "884 884M02 884MF17 884MM18 2 -9\n\
                          884 884MF17 0 0 1 -9\n\
                          884 884MM18 0 0 2 -9\n";

    #[test]
    fn parses_and_counts() {
        let fam = Fam::parse(SAMPLE).unwrap();
        assert_eq!(fam.n_samples(), 3);
        assert_eq!(fam.records[0].iid, "884M02");
        assert_eq!(fam.records[0].sex, "2");
        assert_eq!(fam.records[0].pheno, "-9");
    }

    #[test]
    fn round_trips() {
        let fam = Fam::parse(SAMPLE).unwrap();
        assert_eq!(fam.to_text(), SAMPLE);
    }

    #[test]
    fn zeroes_parents_like_makegeneticmap() {
        let fam = Fam::parse(SAMPLE).unwrap();
        let z = fam.records[0].with_parents_zeroed();
        assert_eq!(z.to_string(), "884 884M02 0 0 2 -9");
    }
}
