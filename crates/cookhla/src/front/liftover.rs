//! Minimal UCSC-chain liftover for chr6 (hg19 → hg18), replacing `pyliftover` for CookHLA's
//! `LiftDown_hg18`. The chr6→chr6 chains are vendored (`data/hg19ToHg18.chr6.chain`) and embedded.
//!
//! CookHLA passes the 1-based `.bim` base position straight to `pyliftover.convert_coordinate`
//! and writes the result back as the new base position; we reproduce that exact behaviour (the
//! input position is used directly as the chain's `t` coordinate). Positions that fall in a chain
//! gap return `None` and the marker is dropped, matching the reference.

/// The vendored chr6→chr6 chains (hg19 target → hg18 query).
const CHR6_CHAIN: &str = include_str!("../../data/hg19ToHg18.chr6.chain");

/// One ungapped aligned block: target `[t, t+size)` maps to query `[q, q+size)`.
#[derive(Debug, Clone, Copy)]
struct Block {
    t: i64,
    q: i64,
    size: i64,
}

#[derive(Debug, Clone)]
struct Chain {
    t_start: i64,
    t_end: i64,
    q_size: i64,
    q_minus: bool,
    blocks: Vec<Block>,
}

/// A chr6 hg19→hg18 lifter.
#[derive(Debug, Clone)]
pub struct LiftOver {
    chains: Vec<Chain>,
}

impl Default for LiftOver {
    fn default() -> Self {
        Self::new()
    }
}

impl LiftOver {
    /// Build from the embedded chr6 chains.
    pub fn new() -> Self {
        Self::parse(CHR6_CHAIN)
    }

    fn parse(text: &str) -> Self {
        let mut chains = Vec::new();
        let mut lines = text.lines().peekable();
        while let Some(line) = lines.next() {
            if !line.starts_with("chain ") {
                continue;
            }
            // chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() < 12 {
                continue;
            }
            let t_start: i64 = f[5].parse().unwrap_or(0);
            let t_end: i64 = f[6].parse().unwrap_or(0);
            let q_size: i64 = f[8].parse().unwrap_or(0);
            let q_minus = f[9] == "-";
            let q_start: i64 = f[10].parse().unwrap_or(0);

            let mut blocks = Vec::new();
            let mut t = t_start;
            let mut q = q_start;
            for bl in lines.by_ref() {
                let b = bl.trim();
                if b.is_empty() {
                    break;
                }
                let parts: Vec<&str> = b.split_whitespace().collect();
                let size: i64 = parts[0].parse().unwrap_or(0);
                blocks.push(Block { t, q, size });
                if parts.len() >= 3 {
                    let dt: i64 = parts[1].parse().unwrap_or(0);
                    let dq: i64 = parts[2].parse().unwrap_or(0);
                    t += size + dt;
                    q += size + dq;
                } else {
                    break; // final block (size only)
                }
            }
            chains.push(Chain {
                t_start,
                t_end,
                q_size,
                q_minus,
                blocks,
            });
        }
        // Prefer larger chains first (they win ties, like pyliftover's score order).
        chains.sort_by_key(|c| -(c.t_end - c.t_start));
        LiftOver { chains }
    }

    /// Lift a chr6 hg19 position to hg18, or `None` if it falls in a gap.
    pub fn convert(&self, pos: i64) -> Option<i64> {
        for c in &self.chains {
            if pos < c.t_start || pos >= c.t_end {
                continue;
            }
            // Binary search the blocks (sorted by t).
            let blocks = &c.blocks;
            let idx = blocks.partition_point(|b| b.t <= pos);
            if idx == 0 {
                continue;
            }
            let b = blocks[idx - 1];
            if pos < b.t + b.size {
                let q = b.q + (pos - b.t);
                return Some(if c.q_minus { c.q_size - 1 - q } else { q });
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lifts_mhc_position_to_reference_coordinate() {
        // Calibrated against the golden: hg19 chr6:29341007 → hg18 29448986 (the reference pos).
        let lo = LiftOver::new();
        assert_eq!(lo.convert(29_341_007), Some(29_448_986));
        // A few more from the example target/reference pairing.
        assert_eq!(lo.convert(29_342_236), Some(29_450_215)); // rs9380122
        assert_eq!(lo.convert(29_342_775), Some(29_450_754)); // rs3749971
    }
}
