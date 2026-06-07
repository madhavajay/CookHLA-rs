#!/usr/bin/env bash
# Drive beagle-rs (the Rust Beagle 5.5 port) for the nine exon×overlap imputations on the real
# CONVERT_IN inputs captured in fixtures/example/golden, writing decompressed VCFs named like the
# oracle's so the consensus parity test can consume them.
#
# Usage:  scripts/run_beaglers_9.sh [out_dir]   (default: fixtures/example/rs)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
G="$ROOT/fixtures/example/golden"
OUT="${1:-$ROOT/fixtures/example/rs}"
BIN="${BEAGLE_RS:-$ROOT/repos/beagle-rs/target/release/beagle-rs}"
ERR="$(head -1 "$ROOT/repos/CookHLA/example/AGM.1958BC+HM_CEU_REF.aver.erate" | tr -d ' ')"

[[ -x "$BIN" ]] || { echo "build beagle-rs first: (cd repos/beagle-rs && cargo build --release -p beagle-rs-cli)" >&2; exit 1; }
[[ -f "$G/1958BC+HM_CEU_REF.MHC.QC.vcf" ]] || { echo "golden CONVERT_IN inputs missing; run docker oracle first" >&2; exit 1; }
mkdir -p "$OUT"

GT="$G/1958BC+HM_CEU_REF.MHC.QC.vcf"
for exon in 2 3 4; do
  REF="$G/HM_CEU_REF.exon${exon}.phased.vcf"
  MAP="$G/AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB.exon${exon}.txt"
  for ol in 0.5 1 1.5; do
    pfx="$OUT/1958BC+HM_CEU_REF.MHC.QC.exon${exon}.${ol}.raw_imputation_out"
    echo ">> beagle-rs exon${exon} overlap ${ol}"
    "$BIN" gt="$GT" ref="$REF" out="$pfx" impute=true gp=true overlap="$ol" \
        err="$ERR" map="$MAP" window=5 ne=10000 nthreads=1 seed=99999 > "$pfx.log" 2>&1
    gzip -dcf "$pfx.vcf.gz" > "$pfx.vcf"
  done
done
echo "beagle-rs imputations in: $OUT"
