//! `beagle2linkage` — convert a Beagle v3 `.bgl` back to a PLINK linkage `.ped` (+ a `.dat`
//! listing the row tags). Native port of `dependency/beagle2linkage.jar` (v1.0).
//!
//! CookHLA pipes a per-exon `.bgl` into it in `HLA_MultipleRefs.Make_ExonN_Panel`
//! (`cat <bgl> | beagle2linkage <prefix>`), then keeps `<prefix>.ped` and discards
//! `<prefix>.dat`. The `.ped` is the transpose of [`super::linkage2beagle`]: one row per
//! sample, the leading columns are the non-`M` header rows (one value per sample), followed by
//! two allele tokens per marker.

use anyhow::{bail, Result};

use crate::bgl::{Bgl, BglLine};

/// Result of `beagle2linkage`: the linkage `.ped` text and the `.dat` row-tag listing.
pub struct LinkageOut {
    pub ped: String,
    pub dat: String,
}

/// Convert a parsed `.bgl` to linkage `.ped` + `.dat` text.
pub fn beagle2linkage(bgl: &Bgl) -> Result<LinkageOut> {
    // Header rows (non-`M`): keep the values after the (tag, label) pair — one per haplotype.
    let headers: Vec<Vec<&str>> = bgl
        .lines
        .iter()
        .filter_map(|l| match l {
            BglLine::Header(raw) => Some(raw.split_ascii_whitespace().skip(2).collect()),
            BglLine::Marker { .. } => None,
        })
        .collect();

    // Marker allele rows.
    let markers: Vec<&Vec<String>> = bgl
        .lines
        .iter()
        .filter_map(|l| match l {
            BglLine::Marker { alleles, .. } => Some(alleles),
            BglLine::Header(_) => None,
        })
        .collect();

    // Haplotype-column count (2 per sample) comes from the first row that has columns.
    let n_hap = headers
        .iter()
        .map(|h| h.len())
        .chain(markers.iter().map(|m| m.len()))
        .next()
        .unwrap_or(0);
    if n_hap % 2 != 0 {
        bail!("beagle2linkage: odd haplotype-column count ({n_hap}); expected diploid pairs");
    }
    let n_samples = n_hap / 2;

    // .ped — one row per sample.
    let mut ped = String::new();
    for s in 0..n_samples {
        let mut first = true;
        let mut push = |val: &str, ped: &mut String| {
            if !first {
                ped.push(' ');
            }
            first = false;
            ped.push_str(val);
        };
        // Leading columns: one value per header row (the first of the sample's haplotype pair).
        for h in &headers {
            push(h.get(2 * s).copied().unwrap_or("0"), &mut ped);
        }
        // Two allele tokens per marker.
        for m in &markers {
            push(m.get(2 * s).map(String::as_str).unwrap_or("0"), &mut ped);
            push(
                m.get(2 * s + 1).map(String::as_str).unwrap_or("0"),
                &mut ped,
            );
        }
        ped.push('\n');
    }

    // .dat — `<tag> \t<label>` for every row (CookHLA discards this, but match the jar anyway).
    let mut dat = String::new();
    for l in &bgl.lines {
        let (tag, label) = match l {
            BglLine::Header(raw) => {
                let mut t = raw.split_ascii_whitespace();
                (t.next().unwrap_or(""), t.next().unwrap_or(""))
            }
            BglLine::Marker { label, .. } => ("M", label.as_str()),
        };
        dat.push_str(tag);
        dat.push_str(" \t"); // jar emits a space then a tab
        dat.push_str(label);
        dat.push('\n');
    }

    Ok(LinkageOut { ped, dat })
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    // Golden captured from `cat l.bgl | java -jar dependency/beagle2linkage.jar bl`.
    const BGL: &str = "P pedigree fam1 fam1 fam2 fam2\n\
                       I id ind1 ind1 ind2 ind2\n\
                       fID father 0 0 0 0\n\
                       mID mother 0 0 0 0\n\
                       C gender 1 1 2 2\n\
                       M rs1 A A A G\n\
                       M rs2 C T T T\n";
    const GOLDEN_PED: &str = "fam1 ind1 0 0 1 A A C T\nfam2 ind2 0 0 2 A G T T\n";
    const GOLDEN_DAT: &str =
        "P \tpedigree\nI \tid\nfID \tfather\nmID \tmother\nC \tgender\nM \trs1\nM \trs2\n";

    #[test]
    fn matches_jar_golden() {
        let bgl = Bgl::parse(BGL).unwrap();
        let out = beagle2linkage(&bgl).unwrap();
        assert_eq!(out.ped, GOLDEN_PED);
        assert_eq!(out.dat, GOLDEN_DAT);
    }

    #[test]
    fn inverts_linkage2beagle() {
        // linkage2beagle(ped) -> bgl ; beagle2linkage(bgl) -> ped should recover the input.
        let ped_in = "fam1 ind1 0 0 1 A A C T\nfam2 ind2 0 0 2 A G T T\n";
        let markers = super::super::linkage2beagle::parse_dat("M rs1\nM rs2\n");
        let bgl_text = super::super::linkage2beagle::linkage2beagle(ped_in, &markers).unwrap();
        let bgl = Bgl::parse(&bgl_text).unwrap();
        assert_eq!(beagle2linkage(&bgl).unwrap().ped, ped_in);
    }
}
