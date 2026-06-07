//! `linkage2beagle` — convert a (no-pheno) PLINK `.ped` + `.dat` marker list into a Beagle v3
//! `.bgl`. Native port of `dependency/linkage2beagle.jar` (v2.0, `standard=true`).
//!
//! CookHLA invokes it in `CONVERT_IN` as
//! `linkage2beagle pedigree=<ped> data=<dat> beagle=<bgl> standard=true`. The `.ped` is one row
//! per sample (`FID IID PID MID SEX` then two allele tokens per marker); the `.dat` lists the
//! markers as `M <id>` rows. The output transposes this to one row per marker, plus the five
//! pedigree header rows Beagle expects.

use anyhow::{bail, Context, Result};

/// Parse a `.dat` (one `M <marker-id>` per line) into the marker ids in order.
pub fn parse_dat(text: &str) -> Vec<String> {
    text.lines()
        .filter_map(|l| {
            let mut t = l.split_ascii_whitespace();
            match (t.next(), t.next()) {
                (Some("M"), Some(id)) => Some(id.to_owned()),
                _ => None,
            }
        })
        .collect()
}

/// Convert `.ped` text + the ordered marker ids to Beagle `.bgl` text.
pub fn linkage2beagle(ped_text: &str, markers: &[String]) -> Result<String> {
    let n_markers = markers.len();

    // Parse each sample row: 5 leading columns then 2 alleles per marker.
    struct Row {
        fid: String,
        iid: String,
        pid: String,
        mid: String,
        sex: String,
        alleles: Vec<String>, // length 2 * n_markers
    }
    let mut rows = Vec::new();
    for (i, line) in ped_text.lines().enumerate() {
        let toks: Vec<&str> = line.split_ascii_whitespace().collect();
        if toks.is_empty() {
            continue;
        }
        if toks.len() < 5 {
            bail!("linkage2beagle: ped row {} has <5 leading columns", i + 1);
        }
        let alleles: Vec<String> = toks[5..].iter().map(|s| s.to_string()).collect();
        if alleles.len() != 2 * n_markers {
            bail!(
                "linkage2beagle: ped row {} has {} allele tokens, expected {} (2 × {} markers)",
                i + 1,
                alleles.len(),
                2 * n_markers,
                n_markers
            );
        }
        rows.push(Row {
            fid: toks[0].into(),
            iid: toks[1].into(),
            pid: toks[2].into(),
            mid: toks[3].into(),
            sex: toks[4].into(),
            alleles,
        });
    }

    let mut out = String::new();
    // Five pedigree header rows, each value repeated once per haplotype column (2× per sample).
    push_header(&mut out, "P pedigree", rows.iter().map(|r| r.fid.as_str()));
    push_header(&mut out, "I id", rows.iter().map(|r| r.iid.as_str()));
    push_header(&mut out, "fID father", rows.iter().map(|r| r.pid.as_str()));
    push_header(&mut out, "mID mother", rows.iter().map(|r| r.mid.as_str()));
    push_header(&mut out, "C gender", rows.iter().map(|r| r.sex.as_str()));

    // One marker row per `.dat` marker, alleles taken column-wise across samples.
    for (j, marker) in markers.iter().enumerate() {
        out.push('M');
        out.push(' ');
        out.push_str(marker);
        for r in &rows {
            let a1 = r
                .alleles
                .get(2 * j)
                .context("linkage2beagle: allele index")?;
            let a2 = r
                .alleles
                .get(2 * j + 1)
                .context("linkage2beagle: allele index")?;
            out.push(' ');
            out.push_str(a1);
            out.push(' ');
            out.push_str(a2);
        }
        out.push('\n');
    }

    Ok(out)
}

/// Write a header row `<tag> <v v>` per sample (each value emitted twice — once per haplotype).
fn push_header<'a>(out: &mut String, tag: &str, vals: impl Iterator<Item = &'a str>) {
    out.push_str(tag);
    for v in vals {
        out.push(' ');
        out.push_str(v);
        out.push(' ');
        out.push_str(v);
    }
    out.push('\n');
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    // Golden captured from
    //   java -jar dependency/linkage2beagle.jar pedigree=l.ped data=l.dat beagle=l.bgl standard=true
    const PED: &str = "fam1 ind1 0 0 1 A A C T\nfam2 ind2 0 0 2 A G T T\n";
    const DAT: &str = "M rs1\nM rs2\n";
    const GOLDEN: &str = "P pedigree fam1 fam1 fam2 fam2\n\
                          I id ind1 ind1 ind2 ind2\n\
                          fID father 0 0 0 0\n\
                          mID mother 0 0 0 0\n\
                          C gender 1 1 2 2\n\
                          M rs1 A A A G\n\
                          M rs2 C T T T\n";

    #[test]
    fn matches_jar_golden() {
        let markers = parse_dat(DAT);
        assert_eq!(markers, ["rs1", "rs2"]);
        assert_eq!(linkage2beagle(PED, &markers).unwrap(), GOLDEN);
    }

    #[test]
    fn round_trips_through_bgl_parser() {
        // The output must parse as a valid Bgl with the right sample ids and markers.
        let markers = parse_dat(DAT);
        let text = linkage2beagle(PED, &markers).unwrap();
        let bgl = crate::bgl::Bgl::parse(&text).unwrap();
        assert_eq!(bgl.sample_ids(), ["ind1", "ind2"]);
        assert_eq!(bgl.marker_labels().collect::<Vec<_>>(), ["rs1", "rs2"]);
    }

    #[test]
    fn rejects_wrong_allele_count() {
        let markers = parse_dat(DAT); // 2 markers => expects 4 allele tokens
        let bad = "fam1 ind1 0 0 1 A A C\n"; // only 3
        assert!(linkage2beagle(bad, &markers).is_err());
    }
}
