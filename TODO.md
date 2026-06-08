# CookHLA → Rust Port — TODO

A Rust reimplementation of [CookHLA](https://github.com/WansonChoi/CookHLA) (accurate HLA
imputation; Cook et al., *Nat Commun* 2021). The original Python/R/csh/Perl source is the
`repos/CookHLA` submodule and is the reference implementation for behavior and tests.

## Status (updated 2026-06-07) — ✅ DEFAULT PIPELINE COMPLETE & VERIFIED END-TO-END

**`cookhla -i 1958BC.hg19 -hg 19 -o OUT -ref HM_CEU_REF -gm <clpsB> -ae <aver.erate>` runs the
complete pipeline from raw input → `.alleles` in ~1.7 s, with HLA calls IDENTICAL to the
original CookHLA.** Every stage is ported to Rust and verified against the captured golden:

| Stage | Module | Verified against golden |
|---|---|---|
| FixInput (liftover hg19→18 + subset + relabel + de-ambiguate) | `front/fixinput.rs` + `front/liftover.rs` (vendored chr6 chain) | raw→NoAmbig, then end-to-end calls |
| QC (MHC extract, strand flip, freq QC, recode) | `front/qc.rs` (+ `plink`) | 1405 markers + reference positions exact |
| Panel building (exon234 + per-exon phased VCFs) | `front/panel.rs` | exon2/3/4 VCFs byte-match (1644/1644/1599) |
| AGM exon-splitting | `front/agm.rs` | end-to-end calls |
| CONVERT_IN (target VCF) | `front/convert_in.rs` | all 1404 `MHC.QC.vcf` rows match |
| Converters (`beagle2vcf`/`linkage2beagle`/`beagle2linkage`) | `convert/` | byte-for-byte vs the jars |
| IMPUTE (`beagle-rs` ×9, parallel) | `impute.rs` | 60 calls, 0 mismatches |
| CONSENSUS (`9GP_no_CI.R`) | `consensus.rs` | 60 calls, 0 mismatches, 1e-15 drift |

**Final parity (`tests/pipeline_parity.rs::raw_end_to_end...`): from raw `1958BC.hg19` →
60/60 HLA calls identical to the oracle golden, 0 mismatches.** Posterior probabilities drift
~few % only because `beagle-rs` is Beagle 5.5 vs the oracle's 5.1; the *calls* match.

**Speed:** ~1.7 s vs the original (Beagle5 imputation alone ≈ 26 s serial + QC/panels/R
consensus/Python orchestration ≈ minutes). No JVM (Beagle runs in-process-class native binary),
9 imputations parallelized with `rayon`, interpreted glue replaced by compiled Rust.

**~4,500 lines of Rust · 50 tests green · `clippy -D warnings` + `fmt` clean.**

- ✅ Workspace + containers + CI. **Multi-arch (amd64 + arm64) `cookhla-rs` Docker image**
  (`docker/Dockerfile.cookhla-rs`): self-contained — `cookhla` + `beagle-rs` (built from the
  submodule) + `plink` (bioconda 1.90b7.7, both arches). **Both arches built and verified to
  produce calls identical to the golden.** CI (`.github/workflows/ci.yml`): `rust`
  (fmt/clippy/test) + `beagle-rs` build + `docker` multi-arch buildx → `ghcr.io/<owner>/cookhla-rs`
  on push/tags. Committed + pushed to `origin/main`.
- ✅ Oracle captured (Beagle5 path): golden `.alleles`/`.hped` + ~200 intermediate stage files.

### ✅ All bundled panels run (adaptive maps) — DONE (2026-06-08)

`cookhla-rs` requires a precomputed adaptive genetic map (`-gm`/`-ae`); the repo only shipped one
for `HM_CEU_REF`, so the 1000 Genomes panels couldn't run. **Fixed by precomputing + bundling a map
per panel** (the `precompute-the-maps` plan — keeps the image multi-arch, no MACH at runtime):
- [x] Ran `MakeGeneticMap` (oracle, deterministic `mach1 seed=123456`, **~24 min/panel**, run 6-way
      in parallel ≈ 25 min; small-sample mode = 200 reference samples → **population-derived, reusable
      maps**) for **all 6 panels** → `.mach_step.avg.clpsB` + `.aver.erate`.
- [x] Bundled in-repo at `data/maps/1000G_REF/<ref_basename>.{mach_step.avg.clpsB,aver.erate}` (~2 MB)
      and in the image at `/opt/cookhla/1000G_REF/maps` (`ENV COOKHLA_MAPS`).
- [x] `cookhla-cli` **auto-resolves** the map by `-ref` basename when `-gm`/`-ae` are omitted
      (`resolve_bundled_map`: `$COOKHLA_MAPS` → image default → local data dir → next to the panel).
- [x] **Verified from a fresh `ghcr.io/.../cookhla-rs:latest` pull: all 6 panels (ALL/AFR/AMR/EAS/
      EUR/SAS) impute with no map args** — 50 calls each (A/B/C/DQB1/DRB1; these panels lack DQA1/DP).
      ALL ~20 s (2504 samples), super-pops ~3–5 s.

Context — how HLA imputation works here: the **panel** (`1000G_REF`) holds reference people's SNPs
*and* HLA alleles; imputation predicts a target's HLA alleles from their SNPs. The **adaptive map**
is a separate per-panel LD-tuning input — a property of the *panel/population*, not the target — so
generating it once per panel and reusing it is correct (and exact for small targets, a strong
approximation for large cohorts).

Still open here:
- [ ] Native `MakeGeneticMap` port (for **non-1000G** references that have no bundled map): STEP0
      randomize + `STEP4-buildMap.R` (`gmap += −½·log(1−rec)`) + `STEP5-collapseHLA.R` + MACH input
      prep. `mach1` is **x86-only** (not in bioconda) → native map-gen would be amd64-only.
- Note: the NYGC high-coverage 1KGP VCFs (GRCh38, per-chromosome) are **not** needed for the above —
  only chr6 matters for HLA and they lack HLA calls (a *fresh* panel would also need HLA typing).

### Other not-yet-ported (later milestones)
- **Legacy paths:** Beagle 4.1 (`-bgl4`), prephasing, `measureAcc`/accuracy, hg38→hg19→hg18
  second liftover hop, `--save-Ambiguous-SNP`, the rare 0-allele `UpdateInput` fix-up.
- **Optimizations:** call `beagle-rs` in-process (lib) instead of subprocess; native PLINK.

Phase checkboxes below reflect this. The remaining work is the milestones above.

## Goal

Make CookHLA **fast and accurate**, then progressively self-contained. CookHLA today is a
thin Python orchestrator that spawns hundreds of subprocesses (PLINK, Beagle, MACH, five
Beagle utility JARs, R, Perl, csh, and awk/sed/grep/cut glue). The compiled binaries it
calls are already fast; the **interpreted orchestration, the R consensus pass, and the
subprocess storms are the slow part**. So:

1. **Keep the fast compiled tools** (PLINK, MACH) as external binaries to start.
2. **Replace the slow interpreted layers in Rust first** — orchestration, the R consensus
   (`9GP_no_CI.R`, an interpreted nested loop over every VCF cell), Perl table-merge, the
   csh driver, and all awk/sed/grep/cut glue.
3. **Run Beagle in-process** via [`beagle-rs`](repos/beagle-rs) — a completed, byte-exact
   Rust port of Beagle 5.5 — instead of shelling out to the JVM nine times.
4. **Replace remaining binaries (PLINK, MACH) natively, incrementally**, once the pipeline
   is green end-to-end.

## Definition of Done (v1 — the fast/accurate milestone)

The default CookHLA path (`Beagle5 + Multiple Markers + Adaptive Genetic Map, no
prephasing`) runs end-to-end in Rust and:

- [ ] Produces `.alleles` + `.hped` with **identical HLA allele calls** to the Python
      reference on the `example/` fixture, with posterior probabilities / confidence within
      a small epsilon (`< 1e-4`). Calls must match exactly; FP noise from MACH and
      accumulation order is tolerated only in the probability fields.
- [ ] Beagle imputation runs **in-process through `beagle-rs`** (no JVM, no nine-fold JVM
      startup), with the nine exon×overlap runs parallelized natively.
- [ ] The Python/R/Perl/csh/awk/sed glue is **gone** — replaced by Rust. PLINK and MACH may
      still be invoked as external binaries (staged replacement, see Phase 11+).
- [ ] **Measurably faster** than the Python baseline on the `example/` fixture (capture a
      before/after wall-clock benchmark; target the orchestration + consensus overhead,
      which is where the interpreted time goes).
- [ ] CLI is **argv-compatible** with `CookHLA.py` for the supported path (`-i -hg -o -ref
      -gm -ae -mem -mp -nth -o`), with matching defaults and validation errors.
- [ ] No `todo!`/`unimplemented!`/panicking stubs in shipped code paths.

**Later milestones** (tracked, not required for v1): port the `MakeGeneticMap`/MACH path,
port the PLINK subset natively, drop all external binaries, support the legacy Beagle4 and
prephasing paths, and the `measureAcc` accuracy module.

---

## Pipeline being ported (the reference data flow)

Default path = `repos/CookHLA/CookHLA.py` → `src/HLA_Imputation_BEAGLE5.py`. Stages:

1. **Input prep** — `src/checkInput.py::FixInput` → liftover target to hg18 (`pyliftover`),
   subset target SNPs to reference by base-position, relabel/flip alleles, drop ambiguous
   `{A,T}`/`{G,C}` SNPs (`exclude_Ambiguous_SNP`). Uses PLINK + pandas + pyliftover.
2. **QC** — `CookHLA.py` (the `FLIP`/`EXTRACT_MHC` blocks): extract MHC `chr6:29–34Mb`,
   strand-flip, allele-freq filter, remove A/T–C/G + non-ACGT, `--geno 0.2`, recode to
   PED/MAP. Uses PLINK + `merge_tables.pl` + awk/sed/cut/grep.
3. **Reference panel + local embedding** — `src/HLA_MultipleRefs.py`,
   `Make_EXON234_Panel.py`, `Make_EXON234_AGM.py`: build the exon234 panel and the per-exon
   (exon2/3/4) sub-panels and AGMs that drive CookHLA's local-embedding trick. Uses PLINK +
   `beagle2linkage.jar` + `beagle2vcf.jar`.
4. **Convert-in** — `CONVERT_IN`: PED→bgl (`linkage2beagle.jar`), refine marker positions
   (`excluding_snp_and_refine_target_position-*.R`), GC-trick P/A↔G/C
   (`bgl2GC_trick_bgl.py`), bgl→VCF (`beagle2vcf.jar`) → target VCF + reference phased VCF.
5. **Impute** — `IMPUTE`: for each `exon ∈ {2,3,4} × overlap ∈ {0.5,1,1.5}` (**9 runs**),
   run **Beagle 5.1** with `gp=true err=<aver.erate> map=<AGM>` → 9 posterior-probability
   VCFs. *This is the engine `beagle-rs` replaces.*
6. **Convert-out / consensus** — `CONVERT_OUT` → `9accuracy_no_CI.v2.csh` → `9GP_no_CI.R`:
   split each VCF per HLA gene; per HLA binary marker combine the 9 imputations
   (dosage `= P(AA) + P(AB)/2`, `max` over the 3 overlaps within an exon, mean across
   exons), normalize to a posterior, pick the top-2 alleles → 2-digit + 4-digit calls +
   posterior probs + confidence. Emit `.alleles`, then `.hped` (`ALLELES2HPED.py`).
7. **Accuracy** *(optional)* — `measureAcc/` vs an answer file.

**Adaptive Genetic Map** (`MakeGeneticMap/`) is a *separate* module: randomize fam, subset
~100–200 samples, run `mach1 --rounds 20 --greedy` to learn recombination/error rates,
build + collapse the map (`STEP4`/`STEP5` R). It produces the `-gm`/`-ae` files. **The
`example/` fixture ships these precomputed** (`example/AGM.1958BC+HM_CEU_REF.{aver.erate,
mach_step.avg.clpsB}`), so the core pipeline can reach end-to-end parity **without MACH** —
MACH/MakeGeneticMap is deferred to Phase 9.

---

## External dependencies — keep / replace plan

| Dependency | Role | Plan |
|---|---|---|
| **Beagle 5.1** (jar) | phasing + imputation HMM (the engine) | **Replace now** — `beagle-rs` lib, in-process |
| **PLINK 1.9** (binary) | bed/bim/fam I/O, QC, freq, flip, recode, subset | **Keep** (fast); port the used subset natively in Phase 11 |
| **MACH** (`mach1` binary) | learn adaptive genetic map | **Keep** (fast); port in Phase 12 (Definition-of-Done-later) |
| linkage2beagle / beagle2linkage / beagle2vcf / vcf2beagle / transpose (jars) | format converters | **Replace now** — small, native Rust (Phase 3) |
| `merge_tables.pl` (Perl) | indexed table join | **Replace now** (Phase 4) |
| `9GP_no_CI.R`, `9accuracy_no_CI.v2.csh` | consensus + HLA calling | **Replace now** — the algorithmic core + a real speed win (Phase 5) |
| `Doubling_vcf.R`, `complete_header.R`, `DP_min_selection.R`, `excluding_snp_*.R`, `STEP0/4/5*.R`, `bgl2GC_trick*.R` | transforms / map build | **Replace** as their stage is reached |
| `checkInput.py` (pandas) + `pyliftover` | input QC + hg→hg18 liftover | **Replace now**; liftover needs UCSC chain files (Phase 8) |
| awk / sed / cut / grep glue | text processing | **Replace now** — native Rust |

Already ported, reused as a library dependency: **`beagle-rs`** (`repos/beagle-rs`, byte-for-byte
Beagle 5.5, drop-in). Other `-rs` ports under `/home/linux/dev/biovault-app/main/repos`
(`samtools-rs`, `bcftools-rs`, `htslib-rs`, `noodles`) are available if VCF/BCF IO helpers
are useful, but Beagle defines the exact VCF/bgl encodings here, so prefer `beagle-rs`'s own
model for anything fed back into it.

---

## Architecture

Cargo workspace, thin-CLI + lib (the convention from `beagle-rs`/`kestrel-rs`):

```
CookHLA-rs/
├── repos/
│   ├── CookHLA/                  # Python reference (submodule, read-only)
│   └── beagle-rs/                # Beagle 5.5 Rust port (submodule, dependency)
├── crates/
│   ├── cookhla/                  # library: pipeline stages, formats, consensus
│   │   ├── src/
│   │   │   ├── plink/            # bed/bim/fam read/write + the op subset (shells to plink first)
│   │   │   ├── bgl/              # .bgl.phased / .markers / .FRQ.frq + GC-trick
│   │   │   ├── convert/          # linkage2beagle / beagle2vcf / vcf2beagle / beagle2linkage / transpose
│   │   │   ├── qc.rs             # MHC extract + strand flip + ambiguous-SNP QC
│   │   │   ├── panel/            # exon234 + per-exon panel & AGM construction
│   │   │   ├── impute.rs         # CONVERT_IN / IMPUTE (beagle-rs) / CONVERT_OUT wiring
│   │   │   ├── consensus.rs      # 9GP_no_CI.R + 9accuracy csh, in Rust
│   │   │   ├── hped.rs           # ALLELES2HPED, allele decoding
│   │   │   ├── agm.rs            # MakeGeneticMap orchestration (mach1 external, Phase 9)
│   │   │   └── lib.rs
│   │   └── tests/                # module parity vs Python intermediates
│   └── cookhla-cli/              # bin: argv-compatible with CookHLA.py
├── fixtures/                     # example/ inputs + golden outputs + intermediates
├── .github/workflows/ci.yml
├── Cargo.toml                    # workspace; beagle-rs as path dep into repos/beagle-rs/crates/beagle-rs
├── TODO.md
└── README.md
```

`beagle-rs` is consumed as a **library** (`beagle_rs::...`, the `main_pkg::main` driver and
the `imp`/`phase` engines) — call it in-process per exon×overlap and parallelize the nine
runs with threads/rayon, instead of `java -jar` × 9.

---

## Guiding principles (from the sibling `-rs` ports)

- **Tests are the spec.** Capture golden outputs from the Python reference and mirror them;
  parity against the reference is the ultimate gate.
- **Parity at every stage boundary, not just the end.** The Python pipeline writes a file at
  each step (`*.MHC.QC.*`, `*.phased.vcf`, the 9 `*.imputation_out.vcf`, `.alleles`). Snapshot
  each as an intermediate golden so a Rust module can be verified in isolation before the
  whole pipeline is wired.
- **Honest gates** (samtools-rs lesson): a missing `plink`/`mach1` must **fail loudly** in a
  preflight check, never silently skip a parity test.
- **Keep fast binaries; replace interpreted glue.** Don't reimplement what's already a fast
  compiled tool until the slow interpreted layer is gone and measured.
- **One focused PR per phase**; `cargo fmt --check` + `clippy -D warnings` + `cargo test`
  green before merge.

---

## Phase 0 — Workspace & harness setup  ✅ DONE (2026-06-07)
- [x] `cargo` workspace: `crates/cookhla` (lib) + `crates/cookhla-cli` (bin). `beagle-rs`
      path dep is stubbed in the workspace manifest, wired at Phase 7 (keeps early builds fast).
- [x] `rust-toolchain.toml` (stable + rustfmt + clippy); base deps (`anyhow`/`thiserror`,
      `rayon`, `clap`, `tempfile`, `pretty_assertions`).
- [x] CI (`.github/workflows/ci.yml`): `fmt --check` + `clippy -D warnings` + `cargo test`.
- [ ] **Preflight**: locate `plink` and `mach1` (PATH or `repos/CookHLA/dependency/`); hard
      error with a clear message if absent. Mirror CookHLA's discovery order. *(Phase 4 — when
      the first `plink` shell-out lands.)*
- [x] `.gitmodules` for both submodules (CookHLA reference + beagle-rs dependency).
- [x] **Containers** (`docker/`): oracle (original CookHLA conda env) + cookhla-rs images,
      `oracle-run.sh`, `parity.sh`, `compose.yml`, `README.md`.

## Phase 1 — Golden oracle & parity harness  ← gate before any porting
- [ ] Stand up the Python reference enough to run the `example/` case end-to-end and emit a
      golden `.alleles` + `.hped`. Command (from README):
      `python CookHLA.py -i example/1958BC.hg19 -hg 19 -o <out> -ref example/HM_CEU_REF
      -gm example/AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB
      -ae example/AGM.1958BC+HM_CEU_REF.aver.erate -mem 2g`.
- [ ] Run with `--save-intermediates` (or patch in the keeps) to capture **every stage
      boundary** as an intermediate golden under `fixtures/example/intermediates/`.
- [ ] Parity comparator: `.alleles`/`.hped` → exact allele-call match + per-field epsilon on
      probabilities. Keep a content-diff of intermediates as a localization diagnostic.
- [ ] Record a **baseline wall-clock** for the Python run (the speed target).
- [ ] Determinism check: does the Python path with the provided AGM reproduce identical
      calls across runs? (Beagle is seeded; with provided AGM, MACH is not invoked.)

## Phase 2 — File formats (bottom-up foundation)  🔄 IN PROGRESS
- [x] PLINK `.bim` reader/writer (tab-out, whitespace-tolerant in, `-1` liftover sentinel).
- [x] PLINK `.fam` reader/writer + `with_parents_zeroed` (the MakeGeneticMap "trick fam") +
      `n_samples` (`getSampleNumbers`).
- [x] PLINK `.frq` reader (header-skipping, padded columns) + `maf_f64`.
- [x] Beagle `.markers` reader/writer + `allele_map`.
- [x] Beagle `.bgl.phased` reader/writer (header/marker rows) + **GC trick** (`Bgl2GC`).
- [ ] PLINK `.bed` genotype matrix (read; write deferred — `plink --make-bed` still does writes).
- [ ] `.ped`/`.map`, `.dat`, `.hped`, `.alleles` models.
- [ ] VCF read/write compatible with what `beagle-rs` ingests/emits (reuse its model where
      possible; this is the handoff format).

## Phase 3 — Format converters (replace the 5 JARs)  🔄 IN PROGRESS
- [x] `beagle2vcf` (`.bgl` + `.markers` → VCF) — **byte-for-byte** vs the jar (fixed
      `##filedate`, `0`→`.` missing handling).
- [x] `linkage2beagle` (PED + `.dat` → `.bgl`, `standard=true`) — byte-for-byte vs the jar.
- [x] GC-trick forward: `bgl2GC_trick_bgl.py::Bgl2GC` (P/A → G/C), tested.
- [ ] `beagle2linkage` (`.bgl` → PED), `vcf2beagle` (VCF → `.bgl`), `transpose`.
- [ ] GC-trick inverse: `GC_tricked_bgl2ori_bgl.py::GCtricedBGL2OriginalBGL`.

## Phase 4 — QC + glue (replace Python/awk/sed/grep/Perl — the orchestration speed win)
- [ ] `merge_tables.pl` → Rust indexed join.
- [ ] CookHLA.py QC blocks (`EXTRACT_MHC`, `FLIP`): MHC region extract, strand-flip,
      `--freq`, frequency merge + parse, A/T–C/G + non-ACGT removal, `--geno 0.2`, recode.
      Orchestrate `plink` calls + native text munging (no per-step subprocess storm).
- [ ] `excluding_snp_and_refine_target_position-*.R`, `redefineBPv1BH.py`,
      `Panel_subset.py`, `SubsetBGLPhased.py` → Rust. Small-sample-mode branch handled.

## Phase 5 — Consensus & HLA calling (the CookHLA-specific algorithm + a speed win)
- [ ] Port `9GP_no_CI.R`: the per-marker, per-sample consensus over the 9 imputation VCFs
      (dosage `P(AA)+P(AB)/2`; `max` over overlaps within an exon; mean across exons;
      column normalize; top-2 allele selection with the `pp/2` second-allele rule;
      2-digit + 4-digit decode; confidence). This interpreted nested loop is a prime
      speedup target — verify call-for-call against the R output.
- [ ] Port `9accuracy_no_CI.v2.csh` orchestration (split each VCF into the 8 HLA genes).
- [ ] `ALLELES2HPED.py`, `Double_alleles_decoder.R`, `complete_header.R`,
      `DP_min_selection.R`, `Doubling_vcf.{R,py}` (only the no-prephasing path needs a subset).

## Phase 6 — Reference panel + local embedding
- [ ] `Make_EXON234_Panel.py`, `Make_EXON234_AGM.py`, `HLA_MultipleRefs.py`
      (`Make_ExonN_Panel`/`Make_ExonN_AGM`): build the exon234 panel + per-exon {2,3,4}
      sub-panels and AGMs; GC-trick + `beagle2vcf` each → `.phased.vcf`. Parallelizable.

## Phase 7 — Imputation wiring (beagle-rs)  🔄 BACK HALF DONE & VERIFIED
- [x] `IMPUTE`: `impute.rs` drives `beagle-rs` per `exon×overlap` with `gp=true err=<aver.erate>
      map=<refined AGM> window ne overlap seed`, nine runs **parallel via rayon**. Verified: from
      real CONVERT_IN inputs, **0/60 call mismatches** vs the Beagle5.1 golden, ~1.5 s.
- [x] `CONVERT_OUT`: feeds the 9 VCFs into the Phase-5 consensus → calls (`impute_and_call`).
- [ ] `CONVERT_IN` (no-prephasing path): assemble the target `MHC.QC.vcf` + reference phased
      VCFs natively (currently consumes the oracle's captured CONVERT_IN inputs). Uses the done
      converters + GC trick + a refine-markers step (`excluding_snp_*.R`, still to port).
- [ ] Optimize: call `beagle-rs` **in-process** (lib) instead of subprocess (later; subprocess
      to a native binary is already fast and avoids coupling to beagle-rs's build).

## Phase 8 — Input prep + liftover
- [ ] `FixInput`/`UpdateInput`: BP-based subset of target→reference, label/allele relabel +
      flip logic (port the `dict_Complement` / 0-allele cases exactly).
- [ ] `LiftDown_hg18`: replace `pyliftover` (hg19→hg18, hg38→hg19→hg18). Vendor the UCSC
      chain files; `chr6` only. (For the hg19 `example/` case, one liftover hop.)
- [ ] `exclude_Ambiguous_SNP`. Wire Phases 2–8 into the full default pipeline.

## Phase 9 — Adaptive Genetic Map (`MakeGeneticMap`) — MACH still external
- [ ] Orchestrate `MakeGeneticMap.py` natively: `STEP0` fam randomize, sample subset,
      format conversions, **invoke `mach1 --rounds 20 --greedy`** (external), `STEP4`
      build-map, `STEP5` collapse → `.aver.erate` + `.mach_step.avg.clpsB`. Small-sample
      mode (target<100 → use reference samples) handled.
- [ ] Wire the `-gm`/`-ae`-absent path so CookHLA-rs can generate its own AGM. Note: MACH's
      RNG means probabilities (not calls) may differ slightly — within the v1 epsilon bar.

## Phase 10 — CLI + end-to-end parity & benchmark
- [ ] `cookhla-cli` argv-compatible with `CookHLA.py` (supported path): args, defaults,
      validation messages, multiprocessing flag.
- [ ] **End-to-end parity** on `example/`: identical HLA calls vs the Phase-1 golden.
- [ ] **Speed benchmark** vs the Python baseline; record the win.

## Phase 11+ — Replace remaining binaries (later milestones)
- [ ] Native PLINK subset (drop the `plink` dependency) — port only the operations used.
- [ ] Native MACH (drop `mach1`) — deterministic port unlocks byte-for-byte parity.
- [ ] Legacy paths: Beagle4.1 (`-bgl4`), prephasing, the `measureAcc`/`NomenCleaner`
      accuracy module, `HapMap` map, `--save-Ambiguous-SNP`.

---

## Decisions (resolved)
- **Scope:** hybrid/staged. Keep PLINK + MACH as fast external binaries first; replace the
  interpreted glue + consensus in Rust now; run Beagle via `beagle-rs` in-process. Replace
  PLINK/MACH natively later.
- **Parity bar (v1):** identical HLA *calls*; probabilities/confidence within `< 1e-4`. Not
  byte-for-byte (that would require a deterministic MACH port, deferred).
- **Primary path:** `Beagle5 + Multiple Markers + AGM, no prephasing` — CookHLA's defaults
  (`f_BEAGLE5=True`, `__use_Multiple_Markers=True`, `f_prephasing=False`). Beagle4 and
  prephasing paths are later milestones.
- **Beagle version:** `beagle-rs` (Beagle 5.5). CookHLA pins 5.1; the `gp=/err=/map=/window=
  /overlap=/ne=` flags used are 5.x-compatible. Validate the CONVERT_IN→IMPUTE handoff
  against the pinned `beagle5.jar` once, then trust beagle-rs's byte-exactness.
- **Golden fixture:** `example/` (`1958BC.hg19` + `HM_CEU_REF`) with the provided AGM — needs
  no MACH, so it isolates the core pipeline for first parity.
- **beagle-rs usage:** library dependency, called in-process; not a shelled binary.

## Open questions
- [ ] How byte-stable is the Python reference's `.alleles` across runs with the provided AGM?
      (Pin the epsilon empirically in Phase 1.)
- [ ] Smallest UCSC chain slice needed for `chr6` hg19/hg38→hg18 liftover to vendor in-repo.
- [ ] Does any consensus tie-break in `9GP_no_CI.R` depend on R's `which.max`/row order in a
      way that affects *calls* (not just probs)? Verify during Phase 5.
- [ ] Should `cookhla` reuse `beagle-rs`'s VCF model directly for intermediates, or keep a
      thin local VCF type? (Decide in Phase 2 once the handoff surface is clear.)

## Reference ports — lessons (`/home/linux/dev/biovault-app/main/repos`)
- **`beagle-rs`** — the dependency here; completed byte-exact Beagle 5.5 port. Two-crate
  workspace, parity-gated against the Java jar, `flate2`/zlib for byte-identical BGZIP.
- **`kestrel-rs`** (Java→Rust) — phased plan, "tests are the spec", preserve documented bugs,
  decide preserve-vs-fix in one final commit.
- **`samtools-rs`** (C→Rust) — upstream parity gate is mandatory; **honest gates** (a missing
  external tool must fail loudly, not silently skip).
- **`bcftools-rs`** (C→Rust) — module-at-a-time, one PR per batch, each slice green on both
  Rust + parity gates before merge.
