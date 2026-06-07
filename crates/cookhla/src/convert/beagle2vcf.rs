//! `beagle2vcf` — convert a phased/unphased `.bgl` (+ its `.markers`) to VCF.
//!
//! Native port of `dependency/beagle2vcf.jar` (v1.2). CookHLA invokes it as
//! `java -jar beagle2vcf.jar <chrom> <markers> <bgl> <missing> > out.vcf`, where `<missing>` is
//! the allele code treated as missing (CookHLA passes `0`, since the GC trick maps unknown
//! calls to `0`). The reference is fed line-for-line by `.markers`, so the `M` rows of the bgl
//! and the `.markers` rows must be 1:1 and in the same order.
//!
//! Output is deterministic — the jar hard-codes `##filedate=20120310` — so we can match it
//! byte-for-byte. Genotype coding: REF = `markers.a1` → `0`, ALT = `markers.a2` → `1`, the
//! missing code → `.`; the two haplotype calls for a sample are joined with `/`.

use anyhow::{bail, Result};

use crate::bgl::{Bgl, Markers};

/// The fixed VCF header the jar emits (note the hard-coded `filedate`).
const HEADER: &str = "##fileformat=VCFv4.1\n\
                      ##filedate=20120310\n\
                      ##source=\"beagle2vcf 1.2\"\n\
                      ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";

/// Convert a `.bgl` + `.markers` to VCF text (tab-separated, trailing newline per row).
///
/// `chrom` is the chromosome label written in column 1 (CookHLA passes `6`). `missing` is the
/// allele string treated as missing in the genotype calls (CookHLA passes `0`).
pub fn beagle2vcf(chrom: &str, markers: &Markers, bgl: &Bgl, missing: &str) -> Result<String> {
    let samples = bgl.sample_ids();
    let marker_rows: Vec<(&str, &Vec<String>)> = bgl
        .lines
        .iter()
        .filter_map(|l| match l {
            crate::bgl::BglLine::Marker { label, alleles } => Some((label.as_str(), alleles)),
            crate::bgl::BglLine::Header(_) => None,
        })
        .collect();

    if marker_rows.len() != markers.len() {
        bail!(
            "beagle2vcf: {} marker rows in bgl but {} in markers file",
            marker_rows.len(),
            markers.len()
        );
    }

    let mut out =
        String::with_capacity(HEADER.len() + marker_rows.len() * (samples.len() * 4 + 40));
    out.push_str(HEADER);

    // Column header line: fixed columns then the sample names.
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for s in &samples {
        out.push('\t');
        out.push_str(s);
    }
    out.push('\n');

    for (i, (label, alleles)) in marker_rows.iter().enumerate() {
        let m = &markers.records[i];
        if *label != m.id {
            bail!(
                "beagle2vcf: marker order mismatch at row {}: bgl {:?} vs markers {:?}",
                i + 1,
                label,
                m.id
            );
        }
        if alleles.len() != samples.len() * 2 {
            bail!(
                "beagle2vcf: marker {label} has {} alleles for {} samples (expected {})",
                alleles.len(),
                samples.len(),
                samples.len() * 2
            );
        }

        // chrom pos id ref alt qual filter info format
        out.push_str(chrom);
        out.push('\t');
        out.push_str(&m.bp.to_string());
        out.push('\t');
        out.push_str(&m.id);
        out.push('\t');
        out.push_str(&m.a1); // REF
        out.push('\t');
        out.push_str(&m.a2); // ALT
        out.push_str("\t.\tPASS\t.\tGT");

        for s in 0..samples.len() {
            let h1 = code(&alleles[2 * s], &m.a1, &m.a2, missing);
            let h2 = code(&alleles[2 * s + 1], &m.a1, &m.a2, missing);
            out.push('\t');
            out.push_str(h1);
            out.push('/');
            out.push_str(h2);
        }
        out.push('\n');
    }

    Ok(out)
}

/// Code one haplotype allele: REF→`0`, ALT→`1`, missing→`.`.
fn code<'a>(allele: &str, ref_a: &str, alt_a: &str, missing: &str) -> &'a str {
    if allele == missing {
        "."
    } else if allele == ref_a {
        "0"
    } else if allele == alt_a {
        "1"
    } else {
        // Not REF/ALT/missing — beagle2vcf has no encoding for this; treat as missing.
        "."
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    // Golden captured from `java -jar dependency/beagle2vcf.jar 6 t.markers t.bgl 0`.
    const BGL: &str = "I id sampleA sampleA sampleB sampleB\n\
                       M rs1 G G C C\n\
                       M rs2 C G G C\n";
    const MARKERS: &str = "rs1 100 G C\nrs2 200 G C\n";
    const GOLDEN: &str = "##fileformat=VCFv4.1\n\
        ##filedate=20120310\n\
        ##source=\"beagle2vcf 1.2\"\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA\tsampleB\n\
        6\t100\trs1\tG\tC\t.\tPASS\t.\tGT\t0/0\t1/1\n\
        6\t200\trs2\tG\tC\t.\tPASS\t.\tGT\t1/0\t0/1\n";

    #[test]
    fn matches_jar_golden() {
        let bgl = Bgl::parse(BGL).unwrap();
        let markers = Markers::parse(MARKERS).unwrap();
        assert_eq!(beagle2vcf("6", &markers, &bgl, "0").unwrap(), GOLDEN);
    }

    #[test]
    fn missing_allele_becomes_dot() {
        // From `... 6 m.markers m.bgl 0` with bgl "M rs1 G 0 C C" -> "0/.\t1/1".
        let bgl = Bgl::parse("I id sA sA sB sB\nM rs1 G 0 C C\n").unwrap();
        let markers = Markers::parse("rs1 100 G C\n").unwrap();
        let vcf = beagle2vcf("6", &markers, &bgl, "0").unwrap();
        let last = vcf.lines().last().unwrap();
        assert_eq!(last, "6\t100\trs1\tG\tC\t.\tPASS\t.\tGT\t0/.\t1/1");
    }

    #[test]
    fn sample_ids_extracted_every_other() {
        let bgl = Bgl::parse(BGL).unwrap();
        assert_eq!(bgl.sample_ids(), ["sampleA", "sampleB"]);
    }
}
