//! Beagle `.markers` — one row per marker: `id  bp  a1  a2` (space-separated).
//!
//! Paired line-for-line with the `M` rows of a `.bgl` file. CookHLA builds it from a `.bim`
//! with `awk '{print $2" "$4" "$5" "$6}'` and feeds it to `beagle2vcf`/the GC trick.

use std::collections::HashMap;
use std::fmt;
use std::path::Path;

use anyhow::{Context, Result};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MarkerRecord {
    pub id: String,
    pub bp: i64,
    pub a1: String,
    pub a2: String,
}

impl MarkerRecord {
    pub fn parse(line: &str) -> Result<Option<Self>> {
        let mut f = crate::ws_fields(line);
        let Some(id) = f.next() else {
            return Ok(None);
        };
        let bp = f.next().context("markers: missing bp (col 2)")?;
        let a1 = f.next().context("markers: missing allele 1 (col 3)")?;
        let a2 = f.next().context("markers: missing allele 2 (col 4)")?;
        let bp: i64 = bp
            .parse()
            .with_context(|| format!("markers: bad bp {bp:?} for {id:?}"))?;
        Ok(Some(MarkerRecord {
            id: id.to_owned(),
            bp,
            a1: a1.to_owned(),
            a2: a2.to_owned(),
        }))
    }
}

impl fmt::Display for MarkerRecord {
    fn fmt(&self, w: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(w, "{} {} {} {}", self.id, self.bp, self.a1, self.a2)
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Markers {
    pub records: Vec<MarkerRecord>,
}

impl Markers {
    pub fn parse(text: &str) -> Result<Self> {
        let mut records = Vec::new();
        for (i, line) in text.lines().enumerate() {
            if let Some(r) = MarkerRecord::parse(line)
                .with_context(|| format!("markers: parse error on line {}", i + 1))?
            {
                records.push(r);
            }
        }
        Ok(Markers { records })
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("markers: reading {}", path.display()))?;
        Self::parse(&text).with_context(|| format!("markers: parsing {}", path.display()))
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
            .with_context(|| format!("markers: writing {}", path.display()))
    }

    /// Map marker id → (a1, a2), for the GC trick's allele lookup.
    pub fn allele_map(&self) -> HashMap<&str, (&str, &str)> {
        self.records
            .iter()
            .map(|r| (r.id.as_str(), (r.a1.as_str(), r.a2.as_str())))
            .collect()
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

    // Real rows from repos/CookHLA/example/HM_CEU_REF.markers.
    const SAMPLE: &str = "rs7772982 29449033 G A\n\
                          rs9380122 29450262 C T\n\
                          rs3749971 29450801 A G\n";

    #[test]
    fn parses_and_round_trips() {
        let m = Markers::parse(SAMPLE).unwrap();
        assert_eq!(m.len(), 3);
        assert_eq!(m.records[1].id, "rs9380122");
        assert_eq!(m.records[1].bp, 29_450_262);
        assert_eq!(m.to_text(), SAMPLE);
    }

    #[test]
    fn allele_map_lookup() {
        let m = Markers::parse(SAMPLE).unwrap();
        let map = m.allele_map();
        assert_eq!(map["rs7772982"], ("G", "A"));
    }
}
