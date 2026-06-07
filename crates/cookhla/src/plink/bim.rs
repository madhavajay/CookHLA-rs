//! PLINK `.bim` — the variant map. Six whitespace-separated columns:
//! `chrom  id  cm  bp  a1  a2`.
//!
//! CookHLA reads/writes `.bim` constantly: cutting columns (`cut -f2,4-`), rebuilding it with
//! awk (`print "6\t" $1 "\t0\t" pos ...`), and (in `checkInput.py`) matching target↔reference
//! markers by base position. `bp` is `i64` because the liftover step uses `-1` as a
//! "failed to lift" sentinel before such markers are dropped.

use std::fmt;
use std::path::Path;

use anyhow::{Context, Result};

/// One `.bim` row. Alleles are kept as owned strings: SNPs are single bases (`A`/`C`/`G`/`T`),
/// but CookHLA's HLA binary markers use `P`/`A` (presence/absence), and `0` marks a missing
/// allele — all must round-trip verbatim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BimRecord {
    pub chrom: String,
    pub id: String,
    pub cm: String, // kept as text: CookHLA always emits "0" here and never does cM math on it
    pub bp: i64,
    pub a1: String,
    pub a2: String,
}

impl BimRecord {
    /// Parse one line. Returns `Ok(None)` for a blank line.
    pub fn parse(line: &str) -> Result<Option<Self>> {
        let mut f = crate::ws_fields(line);
        let Some(chrom) = f.next() else {
            return Ok(None);
        };
        let id = f.next().context("bim: missing id (col 2)")?;
        let cm = f.next().context("bim: missing cM (col 3)")?;
        let bp = f.next().context("bim: missing bp (col 4)")?;
        let a1 = f.next().context("bim: missing allele 1 (col 5)")?;
        let a2 = f.next().context("bim: missing allele 2 (col 6)")?;
        let bp: i64 = bp
            .parse()
            .with_context(|| format!("bim: bad bp {bp:?} for marker {id:?}"))?;
        Ok(Some(BimRecord {
            chrom: chrom.to_owned(),
            id: id.to_owned(),
            cm: cm.to_owned(),
            bp,
            a1: a1.to_owned(),
            a2: a2.to_owned(),
        }))
    }
}

impl fmt::Display for BimRecord {
    /// PLINK writes `.bim` tab-separated; we match that so files are byte-comparable.
    fn fmt(&self, w: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.chrom, self.id, self.cm, self.bp, self.a1, self.a2
        )
    }
}

/// A whole `.bim` file, in order.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Bim {
    pub records: Vec<BimRecord>,
}

impl Bim {
    pub fn parse(text: &str) -> Result<Self> {
        let mut records = Vec::new();
        for (i, line) in text.lines().enumerate() {
            if let Some(r) = BimRecord::parse(line)
                .with_context(|| format!("bim: parse error on line {}", i + 1))?
            {
                records.push(r);
            }
        }
        Ok(Bim { records })
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("bim: reading {}", path.display()))?;
        Self::parse(&text).with_context(|| format!("bim: parsing {}", path.display()))
    }

    /// Serialize back to PLINK's tab-separated text (trailing newline per row).
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
            .with_context(|| format!("bim: writing {}", path.display()))
    }

    pub fn len(&self) -> usize {
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

    // Real rows from repos/CookHLA/example/HM_CEU_REF.bim (tab-separated).
    const SAMPLE: &str = "6\trs7772982\t0\t29448986\tG\tA\n\
                          6\trs9380122\t0\t29450215\tC\tT\n\
                          6\trs3749971\t0\t29450754\tA\tG\n";

    #[test]
    fn parses_example_rows() {
        let bim = Bim::parse(SAMPLE).unwrap();
        assert_eq!(bim.len(), 3);
        assert_eq!(
            bim.records[0],
            BimRecord {
                chrom: "6".into(),
                id: "rs7772982".into(),
                cm: "0".into(),
                bp: 29_448_986,
                a1: "G".into(),
                a2: "A".into(),
            }
        );
        assert_eq!(bim.records[2].bp, 29_450_754);
    }

    #[test]
    fn round_trips_byte_for_byte() {
        let bim = Bim::parse(SAMPLE).unwrap();
        assert_eq!(bim.to_text(), SAMPLE);
    }

    #[test]
    fn tolerates_space_separation_and_blank_lines() {
        let txt = "6 rs1 0 100 A C\n\n6 rs2 0 200 P A\n";
        let bim = Bim::parse(txt).unwrap();
        assert_eq!(bim.len(), 2);
        assert_eq!(bim.records[1].a1, "P"); // HLA presence/absence marker survives
        assert_eq!(bim.records[1].bp, 200);
    }

    #[test]
    fn keeps_liftover_failure_sentinel() {
        let bim = Bim::parse("6 rsX 0 -1 A G\n").unwrap();
        assert_eq!(bim.records[0].bp, -1);
    }
}
