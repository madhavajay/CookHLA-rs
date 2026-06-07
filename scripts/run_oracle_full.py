#!/usr/bin/env python
"""Drive the ORIGINAL CookHLA with __save_intermediates=True so every stage-boundary file is
kept, for stage-by-stage parity testing of cookhla-rs. Runs inside the oracle container.

Usage:  python scripts/run_oracle_full.py <out_prefix>
"""
import sys
from CookHLA import CookHLA

out = sys.argv[1] if len(sys.argv) > 1 else "/golden/1958BC+HM_CEU_REF"

CookHLA(
    "example/1958BC.hg19",   # _input
    "19",                     # _hg_input
    out,                      # _out
    "example/HM_CEU_REF",     # _reference
    "18",                     # _hg_reference
    "example/AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB",  # _AdaptiveGeneticMap (-gm)
    "example/AGM.1958BC+HM_CEU_REF.aver.erate",           # _Average_Erate (-ae)
    _java_memory="2g",
    _MultP=1,
    f_BEAGLE5=True,              # match the CLI default (Beagle 5.1) — what beagle-rs replaces
    __save_intermediates=True,   # keep every intermediate file
)
print("DONE:", out + ".MHC.HLA_IMPUTATION_OUT.alleles")
