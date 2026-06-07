//! `merge_tables.pl` and the small whitespace-table helpers CookHLA's QC builds with
//! awk/sed/cut/grep. Working on structured rows (not byte-exact text) is clearer and matches the
//! semantics the downstream `grep -v -w NA` / `awk` steps rely on.

use anyhow::{bail, Result};

/// A whitespace-delimited table with a header row.
#[derive(Debug, Clone, Default)]
pub struct Table {
    pub header: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl Table {
    /// Parse from text: first non-empty line is the header, the rest are rows (whitespace-split).
    pub fn parse(text: &str) -> Self {
        let mut lines = text.lines().filter(|l| !l.trim().is_empty());
        let header = lines
            .next()
            .map(|l| l.split_whitespace().map(str::to_owned).collect())
            .unwrap_or_default();
        let rows = lines
            .map(|l| l.split_whitespace().map(str::to_owned).collect())
            .collect();
        Table { header, rows }
    }

    pub fn col_index(&self, name: &str) -> Option<usize> {
        self.header.iter().position(|h| h == name)
    }

    /// Keep only rows where no field equals `NA` (the `grep -v -w NA` idiom).
    pub fn drop_na_rows(&mut self) {
        self.rows.retain(|r| !r.iter().any(|f| f == "NA"));
    }
}

/// `merge_tables.pl t1 t2 index`: print every row of `t2`, followed by the matching `t1`
/// non-index columns (joined by the shared `index` value), or `NA` where `t1` has no match.
///
/// Mirrors the Perl: the output header is `t2.header` followed by every `t1` header except the
/// index column; each `t2` row is extended with the matched `t1` values (or `NA`). `t1` rows whose
/// index value is `NA` are not indexed (so they never match).
pub fn merge_tables(t1: &Table, t2: &Table, index: &str) -> Result<Table> {
    let i1 = t1
        .col_index(index)
        .ok_or_else(|| anyhow::anyhow!("merge_tables: index {index:?} not in table 1"))?;
    let i2 = t2
        .col_index(index)
        .ok_or_else(|| anyhow::anyhow!("merge_tables: index {index:?} not in table 2"))?;

    // t1 lookup: index value -> the full t1 row.
    let mut lookup: std::collections::HashMap<&str, &Vec<String>> =
        std::collections::HashMap::new();
    for row in &t1.rows {
        if let Some(key) = row.get(i1) {
            if key != "NA" {
                lookup.entry(key.as_str()).or_insert(row);
            }
        }
    }

    // Columns of t1 to append (all except the index column), with their positions.
    let appended: Vec<usize> = (0..t1.header.len()).filter(|&j| j != i1).collect();

    let mut header = t2.header.clone();
    for &j in &appended {
        header.push(t1.header[j].clone());
    }

    let mut rows = Vec::with_capacity(t2.rows.len());
    for r2 in &t2.rows {
        let mut out = r2.clone();
        let matched = r2.get(i2).and_then(|k| lookup.get(k.as_str()));
        for &j in &appended {
            match matched {
                Some(r1) => out.push(r1.get(j).cloned().unwrap_or_else(|| "NA".into())),
                None => out.push("NA".into()),
            }
        }
        rows.push(out);
    }

    if header.is_empty() {
        bail!("merge_tables: empty header");
    }
    Ok(Table { header, rows })
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn merges_on_index_with_na_for_missing() {
        // t1: reference alleles keyed by SNP; t2: target. Join on SNP.
        let t1 = Table::parse("SNP A1R A2R\nrs1 A G\nrs2 C T\n");
        let t2 = Table::parse("SNP POS A1 A2\nrs1 100 A G\nrs3 300 C G\n");
        let m = merge_tables(&t1, &t2, "SNP").unwrap();
        assert_eq!(m.header, ["SNP", "POS", "A1", "A2", "A1R", "A2R"]);
        // rs1 matched; rs3 not in t1 → NA.
        assert_eq!(m.rows[0], ["rs1", "100", "A", "G", "A", "G"]);
        assert_eq!(m.rows[1], ["rs3", "300", "C", "G", "NA", "NA"]);
    }

    #[test]
    fn drop_na_rows_is_grep_v_w_na() {
        let t1 = Table::parse("SNP A1R A2R\nrs1 A G\n");
        let t2 = Table::parse("SNP POS A1 A2\nrs1 100 A G\nrs3 300 C G\n");
        let mut m = merge_tables(&t1, &t2, "SNP").unwrap();
        m.drop_na_rows();
        assert_eq!(m.rows.len(), 1);
        assert_eq!(m.rows[0][0], "rs1");
    }
}
