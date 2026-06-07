//! Beagle `.bgl`/`.bgl.phased` genotype file + the GC trick.
//!
//! Layout: header rows tagged `P`/`I`/`fID`/`mID`/`C` (pedigree/id/parents/gender) carrying
//! one value per haplotype column, then one `M <marker-id> <allele> <allele> ...` row per
//! marker (two alleles per sample in a phased panel). CookHLA forwards header rows verbatim and
//! rewrites only the `M` rows.

use std::fmt;
use std::path::Path;

use anyhow::{bail, Context, Result};

use super::markers::Markers;

/// One row of a `.bgl` file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BglLine {
    /// A non-`M` header row, preserved verbatim (e.g. `P pedigree ...`, `I id ...`).
    Header(String),
    /// An `M` marker row: the marker label and its per-haplotype allele calls.
    Marker { label: String, alleles: Vec<String> },
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Bgl {
    pub lines: Vec<BglLine>,
}

impl Bgl {
    pub fn parse(text: &str) -> Result<Self> {
        let mut lines = Vec::new();
        for line in text.lines() {
            if line.is_empty() {
                continue;
            }
            // Tag is the first whitespace-delimited token.
            let tag = line.split_ascii_whitespace().next().unwrap_or("");
            if tag == "M" {
                let mut f = line.split_ascii_whitespace();
                f.next(); // "M"
                let label = f
                    .next()
                    .context("bgl: marker row missing label")?
                    .to_owned();
                let alleles = f.map(|s| s.to_owned()).collect();
                lines.push(BglLine::Marker { label, alleles });
            } else {
                lines.push(BglLine::Header(line.to_owned()));
            }
        }
        Ok(Bgl { lines })
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("bgl: reading {}", path.display()))?;
        Self::parse(&text).with_context(|| format!("bgl: parsing {}", path.display()))
    }

    pub fn to_text(&self) -> String {
        let mut s = String::new();
        for l in &self.lines {
            match l {
                BglLine::Header(raw) => s.push_str(raw),
                BglLine::Marker { label, alleles } => {
                    s.push_str("M ");
                    s.push_str(label);
                    for a in alleles {
                        s.push(' ');
                        s.push_str(a);
                    }
                }
            }
            s.push('\n');
        }
        s
    }

    pub fn write(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        std::fs::write(path, self.to_text())
            .with_context(|| format!("bgl: writing {}", path.display()))
    }

    /// Marker labels in file order (the `M` rows), used when subsetting per-exon panels.
    pub fn marker_labels(&self) -> impl Iterator<Item = &str> {
        self.lines.iter().filter_map(|l| match l {
            BglLine::Marker { label, .. } => Some(label.as_str()),
            BglLine::Header(_) => None,
        })
    }

    /// Sample IDs from the `I id ...` header row. The row lists two entries per diploid sample
    /// (one per haplotype column), so the names are every second token (indices 0, 2, 4, ...).
    /// Returns an empty vec if there is no `I id` row.
    pub fn sample_ids(&self) -> Vec<String> {
        for l in &self.lines {
            if let BglLine::Header(raw) = l {
                let mut toks = raw.split_ascii_whitespace();
                if toks.next() == Some("I") && toks.next() == Some("id") {
                    let vals: Vec<&str> = toks.collect();
                    return vals.iter().step_by(2).map(|s| s.to_string()).collect();
                }
            }
        }
        Vec::new()
    }

    /// The GC trick: recode every marker's alleles to `G`/`C` so Beagle's utilities accept them.
    ///
    /// Mirrors `src/bgl2GC_trick_bgl.py::Bgl2GC`: using the `.markers` allele table, each call
    /// equal to A1 becomes `G`, each equal to A2 becomes `C`, anything else (missing) becomes
    /// `0`. Returns the recoded `.bgl` plus the recoded `.markers` (alleles set to `G`/`C`).
    pub fn gc_trick(&self, markers: &Markers) -> Result<(Bgl, Markers)> {
        let amap = markers.allele_map();

        let recoded_lines = self
            .lines
            .iter()
            .map(|l| match l {
                BglLine::Header(raw) => Ok(BglLine::Header(raw.clone())),
                BglLine::Marker { label, alleles } => {
                    let (a1, a2) = amap.get(label.as_str()).copied().ok_or_else(|| {
                        anyhow::anyhow!("Marker ID {label} is not in the marker file")
                    })?;
                    let recoded = alleles
                        .iter()
                        .map(|x| {
                            if x == a1 {
                                "G".to_owned()
                            } else if x == a2 {
                                "C".to_owned()
                            } else {
                                "0".to_owned()
                            }
                        })
                        .collect();
                    Ok(BglLine::Marker {
                        label: label.clone(),
                        alleles: recoded,
                    })
                }
            })
            .collect::<Result<Vec<_>>>()?;

        let recoded_markers = Markers {
            records: markers
                .records
                .iter()
                .map(|r| super::markers::MarkerRecord {
                    id: r.id.clone(),
                    bp: r.bp,
                    a1: "G".into(),
                    a2: "C".into(),
                })
                .collect(),
        };

        // Sanity: a GC-tricked bgl must have one M row per markers row.
        let n_m = recoded_lines
            .iter()
            .filter(|l| matches!(l, BglLine::Marker { .. }))
            .count();
        if n_m != recoded_markers.len() {
            bail!(
                "gc_trick: {} marker rows in bgl but {} in markers file",
                n_m,
                recoded_markers.len()
            );
        }

        Ok((
            Bgl {
                lines: recoded_lines,
            },
            recoded_markers,
        ))
    }
}

impl fmt::Display for Bgl {
    fn fmt(&self, w: &mut fmt::Formatter<'_>) -> fmt::Result {
        w.write_str(&self.to_text())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bgl::markers::Markers;
    use pretty_assertions::assert_eq;

    // Minimal phased bgl: 2 samples (4 haplotype columns), 2 markers — same shape as the
    // reference HM_CEU_REF.bgl.phased, just trimmed.
    const BGL: &str = "P pedigree s1 s1 s2 s2\n\
                       I id s1 s1 s2 s2\n\
                       M rs1 A A G G\n\
                       M rs2 T C C T\n";

    const MARKERS: &str = "rs1 100 A G\n\
                           rs2 200 C T\n";

    #[test]
    fn parses_headers_and_markers() {
        let bgl = Bgl::parse(BGL).unwrap();
        assert_eq!(bgl.lines.len(), 4);
        assert!(matches!(bgl.lines[0], BglLine::Header(_)));
        let labels: Vec<_> = bgl.marker_labels().collect();
        assert_eq!(labels, ["rs1", "rs2"]);
    }

    #[test]
    fn round_trips() {
        assert_eq!(Bgl::parse(BGL).unwrap().to_text(), BGL);
    }

    #[test]
    fn gc_trick_recodes_alleles_and_markers() {
        let bgl = Bgl::parse(BGL).unwrap();
        let markers = Markers::parse(MARKERS).unwrap();
        let (gc_bgl, gc_markers) = bgl.gc_trick(&markers).unwrap();

        // rs1: A1=A→G, A2=G→C  => "A A G G" -> "G G C C"
        // rs2: A1=C→G, A2=T→C  => "T C C T" -> "C G G C"
        let expected = "P pedigree s1 s1 s2 s2\n\
                        I id s1 s1 s2 s2\n\
                        M rs1 G G C C\n\
                        M rs2 C G G C\n";
        assert_eq!(gc_bgl.to_text(), expected);

        // Markers recoded to G/C, positions preserved.
        assert_eq!(gc_markers.to_text(), "rs1 100 G C\nrs2 200 G C\n");
    }

    #[test]
    fn gc_trick_maps_unknown_allele_to_zero() {
        // A genotype call that matches neither A1 nor A2 (e.g. missing "0") becomes "0".
        let bgl = Bgl::parse("M rs1 A 0 G N\n").unwrap();
        let markers = Markers::parse("rs1 100 A G\n").unwrap();
        let (gc_bgl, _) = bgl.gc_trick(&markers).unwrap();
        assert_eq!(gc_bgl.to_text(), "M rs1 G 0 C 0\n");
    }
}
