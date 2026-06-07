//! Consensus HLA caller — the heart of `CONVERT_OUT`, ported from `src/9GP_no_CI.R` (driven by
//! `src/9accuracy_no_CI.v2.csh`). This is CookHLA's local-embedding combiner **and** the single
//! biggest speed win: the R original is an interpreted nested loop over every (marker × sample)
//! cell of nine VCFs; here it is plain compiled float math.
//!
//! Input: the nine imputation-output VCFs = `exon ∈ {2,3,4} × overlap (3 each)`, ordered by exon
//! then overlap. For each HLA gene and each binary allele marker (`HLA_<gene>_<allele>_exon<N>`):
//!
//! 1. Per imputation, per sample, the value is `GP[0] + GP[1]/2` (the `GP` FORMAT subfield —
//!    `P(0/0) + P(0/1)/2`, the expected reference-allele dose / 2).
//! 2. Combine the nine: `max` over the three overlaps **within** each exon (best overlap), then
//!    average across exons — `(m2 + m3 + m4)/3`, or `(m2 + m3)/2` for alleles the exon-4 panel
//!    lacks (the R `if i <= nrow(HLA_EXON7)` guard).
//! 3. Normalize each sample's values to a posterior over the gene's alleles.
//! 4. Call the top allele, then the best second allele with posterior `> pp1/2` (else homozygous);
//!    confidence = `pp1` (hom) or `pp1 + pp2` (het). Decode 2-digit + 4-digit from the marker id.

use std::collections::BTreeMap;

use anyhow::Result;

use crate::HLA_NAMES;

/// One imputation VCF reduced to its HLA allele rows + sample order.
#[derive(Debug, Clone)]
pub struct ImputationVcf {
    pub samples: Vec<String>,
    /// HLA allele rows in file order: `(marker_id, per-sample value = GP0 + GP1/2)`.
    pub hla_rows: Vec<(String, Vec<f64>)>,
}

impl ImputationVcf {
    /// Parse a Beagle imputation-output VCF (`GT:DS:GP`), keeping only `HLA_*` rows and reducing
    /// each genotype cell to `GP[0] + GP[1]/2`. Replaces the csh `grep "#CHROM"` + `grep HLA_*`.
    pub fn parse(text: &str) -> Result<Self> {
        let mut samples = Vec::new();
        let mut hla_rows = Vec::new();

        for line in text.lines() {
            if line.starts_with("##") {
                continue;
            }
            if line.starts_with("#CHROM") || line.starts_with("CHROM") {
                // Header: columns 9.. (0-based) are sample names.
                samples = line
                    .split(['\t', ' '])
                    .filter(|s| !s.is_empty())
                    .skip(9)
                    .map(|s| s.to_owned())
                    .collect();
                continue;
            }
            // Data row.
            let mut fields = line.split(['\t', ' ']).filter(|s| !s.is_empty());
            let cols: Vec<&str> = fields.by_ref().take(9).collect();
            if cols.len() < 9 {
                continue;
            }
            let id = cols[2];
            if !id.starts_with("HLA_") {
                continue;
            }
            let gp_idx = gp_index(cols[8]);
            let values: Vec<f64> = fields.map(|cell| cell_value(cell, gp_idx)).collect();
            hla_rows.push((id.to_owned(), values));
        }

        Ok(ImputationVcf { samples, hla_rows })
    }
}

/// Index of the `GP` subfield in a `FORMAT` string (e.g. `GT:DS:GP` → 2). Defaults to 2.
fn gp_index(format: &str) -> usize {
    format.split(':').position(|f| f == "GP").unwrap_or(2)
}

/// `GP[0] + GP[1]/2` from a genotype cell like `0|1:0.87:0.13,0.87,0`.
fn cell_value(cell: &str, gp_idx: usize) -> f64 {
    let Some(gp) = cell.split(':').nth(gp_idx) else {
        return 0.0;
    };
    let mut it = gp.split(',');
    let p0: f64 = it.next().and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let p1: f64 = it.next().and_then(|s| s.parse().ok()).unwrap_or(0.0);
    p0 + p1 / 2.0
}

/// The final call for one (sample, gene).
#[derive(Debug, Clone, PartialEq)]
pub struct AlleleCall {
    pub fid: String,
    pub gene: String,
    pub two_digit: (String, String),
    pub four_digit: (String, String),
    pub pp1: f64,
    pub pp2: f64,
    pub conf: f64,
}

impl AlleleCall {
    /// `.alleles` line: `FID IID gene 2d1,2d2 4d1,4d2 pp1 pp2 conf` (space-separated).
    pub fn to_alleles_line(&self) -> String {
        format!(
            "{fid} {fid} {gene} {a2a},{a2b} {a4a},{a4b} {pp1} {pp2} {conf}",
            fid = self.fid,
            gene = self.gene,
            a2a = self.two_digit.0,
            a2b = self.two_digit.1,
            a4a = self.four_digit.0,
            a4b = self.four_digit.1,
            pp1 = fmt_prob(self.pp1),
            pp2 = fmt_prob(self.pp2),
            conf = fmt_prob(self.conf),
        )
    }
}

/// Format a probability like R's default `write.table` (drop trailing zeros; `0.5`, `1`, not
/// `0.500000`). Calls must match exactly; probabilities only need to match within epsilon, so we
/// render up to 15 significant digits and trim.
fn fmt_prob(x: f64) -> String {
    if x == x.trunc() && x.abs() < 1e15 {
        return format!("{}", x as i64);
    }
    let s = format!("{x:.15}");
    let s = s.trim_end_matches('0');
    s.trim_end_matches('.').to_owned()
}

/// Decode the allele-digit string from a marker id: `HLA_A_0101_exon2` → `0101`.
fn allele_digits(id: &str) -> &str {
    id.split('_').nth(2).unwrap_or("")
}

/// Run the consensus for a single gene over the nine imputation VCFs (ordered exon2×3, exon3×3,
/// exon4×3). Returns one call per sample, or an empty vec if the gene is absent.
// The combine/normalize passes are genuine 2D (allele × sample) grid walks; indexing both axes
// is clearer than nested iterator zips here.
#[allow(clippy::needless_range_loop)]
pub fn call_gene(gene: &str, nine: &[ImputationVcf; 9]) -> Vec<AlleleCall> {
    let prefix = format!("HLA_{gene}_");

    // Per file, the gene's allele rows in order.
    let per_file: Vec<Vec<&(String, Vec<f64>)>> = nine
        .iter()
        .map(|v| {
            v.hla_rows
                .iter()
                .filter(|(id, _)| id.starts_with(&prefix))
                .collect()
        })
        .collect();

    let n_alleles = per_file.iter().map(|r| r.len()).max().unwrap_or(0);
    if n_alleles == 0 {
        return Vec::new(); // gene absent in this study
    }
    let samples = &nine[0].samples;
    let n_samples = samples.len();
    let n_exon4 = per_file[6].len(); // rows present in the exon-4 panel

    // value(file, allele_row, sample), 0 if out of range.
    let val = |f: usize, r: usize, s: usize| -> f64 {
        per_file[f]
            .get(r)
            .and_then(|(_, vals)| vals.get(s).copied())
            .unwrap_or(0.0)
    };

    // Combined (un-normalized) posterior weight per (allele, sample).
    let mut combined = vec![vec![0.0f64; n_samples]; n_alleles];
    for r in 0..n_alleles {
        for s in 0..n_samples {
            let m2 = val(0, r, s).max(val(1, r, s)).max(val(2, r, s));
            let m3 = val(3, r, s).max(val(4, r, s)).max(val(5, r, s));
            combined[r][s] = if r < n_exon4 {
                let m4 = val(6, r, s).max(val(7, r, s)).max(val(8, r, s));
                (m2 + m3 + m4) / 3.0
            } else {
                (m2 + m3) / 2.0
            };
        }
    }

    // Allele ids in row order (from the first file, suffix-independent for decoding).
    let ids: Vec<&str> = (0..n_alleles)
        .map(|r| per_file[0].get(r).map(|(id, _)| id.as_str()).unwrap_or(""))
        .collect();

    // Normalize each sample column to a posterior over alleles.
    for s in 0..n_samples {
        let colsum: f64 = (0..n_alleles).map(|r| combined[r][s]).sum();
        if colsum > 0.0 {
            for r in 0..n_alleles {
                combined[r][s] /= colsum;
            }
        }
    }

    // Call per sample.
    let mut calls = Vec::with_capacity(n_samples);
    for s in 0..n_samples {
        // Top allele (first argmax).
        let mut r1 = 0;
        for r in 1..n_alleles {
            if combined[r][s] > combined[r1][s] {
                r1 = r;
            }
        }
        let pp1 = combined[r1][s];
        let id1 = ids[r1];

        // Best second allele with posterior > pp1/2 (id distinct from the first).
        let mut id2 = id1;
        let mut pp2 = pp1;
        let mut conf = pp1;
        let mut best = 0.0;
        let mut found = false;
        for r in 0..n_alleles {
            let v = combined[r][s];
            if ids[r] != id1 && v > pp1 / 2.0 && v > best {
                id2 = ids[r];
                pp2 = v;
                conf = pp1 + v;
                best = v;
                found = true;
            }
        }
        let _ = found;

        let d1 = allele_digits(id1);
        let d2 = allele_digits(id2);
        calls.push(AlleleCall {
            fid: samples[s].clone(),
            gene: gene.to_owned(),
            two_digit: (two_digit(d1), two_digit(d2)),
            four_digit: (d1.to_owned(), d2.to_owned()),
            pp1,
            pp2,
            conf,
        });
    }
    calls
}

/// First two characters of an allele-digit string (the 1-field / 2-digit allele).
fn two_digit(digits: &str) -> String {
    digits.chars().take(2).collect()
}

/// Full `CONVERT_OUT`: run every HLA gene and concatenate the `.alleles` rows (gene order matches
/// the reference: A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1). Genes absent from the data are skipped.
pub fn convert_out(nine: &[ImputationVcf; 9]) -> Vec<AlleleCall> {
    // Gene output order in CookHLA's csh is A B C DRB1 DQA1 DQB1 DPA1 DPB1; the merge step then
    // re-emits in HLA_NAMES order. We follow HLA_NAMES (A B C DPA1 DPB1 DQA1 DQB1 DRB1).
    let mut out = Vec::new();
    for gene in HLA_NAMES {
        out.extend(call_gene(gene, nine));
    }
    out
}

/// Group nine imputation VCFs keyed by `(exon, overlap)` into the exon-then-overlap order the
/// consensus expects. `exon` ∈ {2,3,4}; overlaps are sorted ascending within each exon (the
/// within-exon `max` makes the order immaterial, but we sort for determinism).
pub fn order_nine(
    by_key: BTreeMap<(u8, String), ImputationVcf>,
) -> Result<[ImputationVcf; 9], anyhow::Error> {
    let mut v: Vec<ImputationVcf> = Vec::with_capacity(9);
    for exon in [2u8, 3, 4] {
        let mut this: Vec<_> = by_key.iter().filter(|((e, _), _)| *e == exon).collect();
        this.sort_by(|a, b| a.0 .1.cmp(&b.0 .1));
        for (_, vcf) in this {
            v.push(vcf.clone());
        }
    }
    if v.len() != 9 {
        anyhow::bail!("order_nine: expected 9 imputation VCFs, got {}", v.len());
    }
    v.try_into()
        .map_err(|_| anyhow::anyhow!("order_nine: could not collect 9 VCFs"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vcf(rows: &[(&str, &[f64])]) -> ImputationVcf {
        ImputationVcf {
            samples: vec!["S1".into()],
            hla_rows: rows
                .iter()
                .map(|(id, vals)| (id.to_string(), vals.to_vec()))
                .collect(),
        }
    }

    #[test]
    fn gp_value_extraction() {
        // GP = field index 2; value = 0.13 + 0.87/2 = 0.565.
        assert!((cell_value("0|1:0.87:0.13,0.87,0", 2) - 0.565).abs() < 1e-12);
    }

    #[test]
    fn parses_hla_rows_only() {
        let text = "##fileformat=VCFv4.2\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
            6\t100\trs1\tG\tC\t.\tPASS\t.\tGT:DS:GP\t0|0:0:1,0,0\n\
            6\t200\tHLA_A_0101_exon2\tG\tC\t.\tPASS\t.\tGT:DS:GP\t0|1:1:0.2,0.8,0\n";
        let v = ImputationVcf::parse(text).unwrap();
        assert_eq!(v.samples, ["S1"]);
        assert_eq!(v.hla_rows.len(), 1); // rs1 skipped
        assert_eq!(v.hla_rows[0].0, "HLA_A_0101_exon2");
        assert!((v.hla_rows[0].1[0] - (0.2 + 0.8 / 2.0)).abs() < 1e-12);
    }

    #[test]
    fn homozygous_call_when_one_allele_dominates() {
        // One allele with value 1 everywhere, another with 0. After normalize: pp1=1, no second.
        let strong = ("HLA_A_0101_exon2", &[1.0f64][..]);
        let weak = ("HLA_A_0201_exon2", &[0.0f64][..]);
        let nine: [ImputationVcf; 9] = std::array::from_fn(|_| vcf(&[strong, weak]));
        let calls = call_gene("A", &nine);
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].four_digit, ("0101".into(), "0101".into()));
        assert_eq!(calls[0].two_digit, ("01".into(), "01".into()));
        assert!((calls[0].pp1 - 1.0).abs() < 1e-12);
        assert!((calls[0].conf - 1.0).abs() < 1e-12); // homozygous: conf = pp1
    }

    #[test]
    fn heterozygous_call_picks_two_alleles() {
        // Two alleles each value 1 → normalized 0.5/0.5 → het, conf = 1.
        let a = ("HLA_A_0101_exon2", &[1.0f64][..]);
        let b = ("HLA_A_0201_exon2", &[1.0f64][..]);
        let nine: [ImputationVcf; 9] = std::array::from_fn(|_| vcf(&[a, b]));
        let calls = call_gene("A", &nine);
        assert_eq!(calls[0].four_digit, ("0101".into(), "0201".into()));
        assert!((calls[0].pp1 - 0.5).abs() < 1e-12);
        assert!((calls[0].pp2 - 0.5).abs() < 1e-12);
        assert!((calls[0].conf - 1.0).abs() < 1e-12);
    }

    #[test]
    fn fmt_prob_matches_r_style() {
        assert_eq!(fmt_prob(1.0), "1");
        assert_eq!(fmt_prob(0.5), "0.5");
        assert_eq!(fmt_prob(0.975), "0.975");
    }
}
