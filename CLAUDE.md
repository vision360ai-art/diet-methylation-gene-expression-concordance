# CLAUDE.md

Context for future sessions. Read this first — it captures everything non-obvious
so you don't have to re-derive it. (Raw data is gitignored, so some facts here are
not visible from the repo itself.)

## Research question
In the LIFE-Adult cohort (Leipzig, blood), which **CpG methylation** changes
associated with **diet quality** (continuous `DietScore`) show **concordant
changes in expression of the corresponding gene** — independent of BMI, age, and
sex? This applies an established methylation↔expression concordance framework to a
dataset/question combination that does not appear to be previously published
(extends Klemp et al. 2022, which analyzed LIFE-Adult methylation only).

> **⚠️ EXPRESSION-SIDE PIVOT (2026-07-11).** The paired LIFE-Adult expression matrix
> never became available. Expression now comes from a **separate public cohort,
> GSE109597** (see Data structure). This breaks the original *within-subject* design:
> methylation and expression are no longer from the same people, so the **eQTM /
> paired concordance stage cannot be computed within samples** as written. The
> question must be reframed as a **cross-cohort, cross-platform, cross-exposure**
> integration (diet→methylation in LIFE vs. obesity→expression in GSE109597), with
> the mismatch documented honestly. Working out that reframing + adapting `02`/`03` is
> the main next-session task — it is NOT yet resolved in the code.

## Data structure (`data/raw/`, gitignored)
| File | What | n | Platform |
|---|---|---|---|
| `SampleSheet_V2.csv` | Methylation sample sheet + phenotypes | 98 | Illumina MethylationEPIC (~850K), blood |
| `Gx_sample_info-21-06-15_v3.txt` | Expression sample metadata | 1,883 | Illumina HumanHT-12 v4, blood |
| `Gx_probe_info-21-06-15_v3.txt` | Expression probe annotation + QC flags (48,106 probes × 42 cols) | — | — |
| `Gx_DataDictionary-21-06-15_v3.xlsx` | Column docs for sample_info + probe annotation | — | — |
| `Clinical…Klemp…2022.pdf` | Source paper (candidate-gene provenance, see below) | — | — |
| `GSE109597_series_matrix.txt.gz` | **Expression VALUES — the actual expression source now** | 84 | **Affymetrix U133 Plus 2.0 (GPL570)**, whole blood |

- **The LIFE-Adult `Gx_*` files are methylation-paired but VALUE-LESS:** only probe
  annotation + sample metadata ever arrived, never the per-sample values. They are now
  effectively legacy — the expression *values* come from **GSE109597** instead.
- **GSE109597** ("Predictive computational obesity risk framework…", NIH/NINR): 84
  whole-blood samples (design says n=90 — reconcile the 6-sample gap next session),
  Affymetrix U133 Plus 2.0, **already log2/RMA-normalized**. Value matrix is `ID_REF`
  (Affy probeset, e.g. `1007_s_at`; 54,675 probesets) × `GSM…` sample columns,
  starting after `!series_matrix_table_begin` (line 71). Exposure framing is
  **obesity/BMI risk**, not diet — different people AND a different question.
- **Cross-assay linkage NO LONGER HOLDS for expression.** The old
  `Synonym`==`sampleID` (48 paired LIFE subjects, 25F/23M, `DietScore` 3–24) linked
  methylation to the LIFE expression that never materialised; GSE109597 samples are
  unrelated individuals, so there is no within-subject key to join on.
- **Tracked curated input:** `data/raw/annotation/candidate_genes.txt` is the one
  file under `data/raw/` that IS git-tracked (a `.gitignore` exception) — it's a
  derived analysis input, not bulk raw data.

## Known data quirks
- **Duplicate `HLS` column** in `SampleSheet_V2.csv`: appears twice — categorical
  (`low`/`high`) and numeric (composite Healthy Lifestyle Score). Load with
  `read.csv(check.names = FALSE)` and rename explicitly to **`HLS_group`** and
  **`HLS_score`** (done in `01_preprocess_methylation.R`).
- **IDAT files not yet in repo.** minfi needs the raw `*_Grn.idat`/`*_Red.idat`
  pairs; they are gitignored and expected in **`data/raw/idat/`**. Script 01 fails
  fast with a clear message if absent.
- **Expression source = GSE109597 (Affymetrix), NOT the LIFE Illumina data.**
  `02_preprocess_expression.R` as written targets the LIFE HT-12 layout (Illumina
  probe IDs, `symbol_INGENUITY` gene column, provider QC flags) — **none of that
  applies to GSE109597**. The GSE matrix is already log2/RMA-normalized (so no
  provider-flag filtering step), keyed by Affy probeset IDs that need a **GPL570 /
  `hgu133plus2.db` probeset→gene-symbol mapping** (a different annotation than
  Illumina), with multi-probeset-per-gene collapse. `02` needs a real rewrite for this
  source — the `# CONFIRM` items are moot. This is next-session work.

## Current status & blockers (2026-07-11)
**Methylation side is RUN and the EWAS is done.** `01` has run on the **36** paired
samples that have usable IDATs (12 of the 48 IDAT pairs are still missing) →
`data/processed/{beta,m}_values.rds`, `grSet…rds`, `pheno_methylation.{rds,csv}` with
the 6 Houseman cell-type proportions. Cell-type deconvolution runs via the low-memory
`01a_cell_counts.R` path (see [[houseman-lowmem-and-host-limits]] memory): the host
OOMs on `estimateCellCounts2`, so 01 reads a cached `cell_counts.csv`.

**EWAS result** (`03a_ewas.R`, the carved-out Sections 0–3 of `03`; model
`M ~ DietScore + Age + Sex + BMI + 5 cell-type props`, n=36, 806,845 CpGs):
- **0** CpGs at genome-wide FDR<0.05 — *expected* at this n (EWAS is discovery, not
  the inference stage).
- **996** CpGs at nominal **p<1e-3** (100 at p<1e-4, 9 at p<1e-5). Top hit
  **cg07018629, p=6.9e-7**.
- These 996 are the **discovery pool** that feeds the concordance analysis (unioned
  with top-N, then FDR-controlled *within* the CpG–gene pair set) once expression
  data arrives. Tracked as `results/tables/ewas_diet_significant.csv`; the full
  806k-row table is gitignored (`ewas_diet_methylation.csv.gz`, regenerable).

**Expression data now IN HAND — but from a different cohort (see the pivot note at
top).** `GSE109597_series_matrix.txt.gz` (Affymetrix U133 Plus 2.0, 84 whole-blood
samples, obesity-risk study) is downloaded + verified. So the remaining work is not
"download" but "**reconcile**": (a) rewrite `02` to parse the GSE series matrix (Affy
probeset→gene via GPL570/`hgu133plus2.db`; already log2/RMA-normalized, so no
provider-flag step); (b) redesign `03`'s concordance since the eQTM cannot be
within-subject across two unrelated cohorts; (c) document the cohort/platform/exposure
mismatch in the methods. Full plain-language status lives in
[`docs/execution_plan.md`](docs/execution_plan.md); this file is the technical brief.

**Exploratory QC** (`explore_methylation_structure.R`, raw vs. cell-corrected): raw
methylation PC1 is almost entirely blood cell composition (PC1↔neutrophil r=−0.95);
residualizing on the 6 proportions removes it and **age** emerges as the next axis
(PC2↔Age r=−0.59). DietScore shows no global PCA/clustering structure before or after
correction — expected, and the motivation for the supervised, covariate-adjusted EWAS.
Figures: `results/figures/pca_{raw,cellcorr}_*`, `hclust_{raw,cellcorr}_*`.

**Next session:** (1) quick + independent — wire `maxprobes` cross-reactive filtering
into `01`. (2) the major task — parse the `GSE109597` value section (past the line-71
`!series_matrix_table_begin` header), rewrite `02` for Affymetrix/GSE109597, and
redesign `03`'s concordance for the cross-cohort / cross-exposure reality, with the
mismatch caveat written into the methods (see [[gse109597-expression-pivot]]).

## Known limitations (stated upfront, not compromises)
- **n=48 paired is underpowered** for genome-wide FDR across ~850K probes (and only
  **36** are currently usable — 12 IDAT pairs still missing; the EWAS ran on 36).
  Handled by treating EWAS as *discovery* (nominal-p + top-N) and controlling FDR
  *within* the CpG–gene pair set, plus a candidate-gene arm.
- **Both platforms are arrays** (methylation EPIC; expression now **Affymetrix U133
  Plus 2.0** via GSE109597, not the originally-planned Illumina HT-12v4), not
  WGBS/RRBS or RNA-seq.
- **Cross-cohort / cross-exposure expression — the biggest caveat, from the pivot.**
  Expression (GSE109597) and methylation (LIFE-Adult) are DIFFERENT people, DIFFERENT
  array platforms, and DIFFERENT exposures (obesity/BMI risk vs. diet quality). There
  is no within-subject eQTM; "concordance" degrades to gene-level agreement in effect
  *direction* across two cohorts — much weaker, and must be stated plainly, not
  papered over. Reframing the analysis honestly around this is next-session work.
- **Blood cell composition** is the dominant methylation confounder — corrected via
  Houseman deconvolution (see pipeline).

## Analysis plan / pipeline (`scripts/`, numbered by order)
1. **`01_preprocess_methylation.R`** — minfi: read IDATs → QC (detP, sex check) →
   drop failed samples → **Houseman cell-type deconvolution** (FlowSorted.Blood.EPIC/
   IDOL; 6 proportions written into pheno) → **noob+quantile normalization**
   (`cfg$normalization` default — funnorm is available as a sensitivity option but has
   the highest peak memory and OOMs on this host) → probe filtering → Beta/M matrices.
   NB: deconvolution is done by **`01a_cell_counts.R`** (low-memory
   `projectCellType_CP` + compTable path) and cached to `cell_counts.csv`; the inline
   `estimateCellCounts2` OOMs on this host and is only a fallback.
   `01b_finalize_from_grset.R` recovers the save step if it fails midway. The
   genome-wide EWAS itself lives in **`03a_ewas.R`** (Sections 0–3 of `03`), runnable
   without expression data.
   **Cross-reactive probe filtering currently SELF-SKIPS** — it looks for
   `data/raw/annotation/Pidsley2016_crossreactive.csv` (absent). The **`maxprobes`**
   package (github.com/markgene/maxprobes; bundles the Pidsley et al. 2016 EPIC
   cross-reactive list, 43,256 probes) is now installed to supply it but is **NOT yet
   wired in — that is the FIRST task next session** (see [[maxprobes-crossreactive]]:
   arg is `array_type`, returns a `list` so `unlist()` before `%in%`; depends on
   `minfiData`).
2. **`02_preprocess_expression.R`** — **STALE / needs rewrite for GSE109597.** As
   written it loads the LIFE HT-12 matrix → provider QC-flag filtering
   (`probe_QCok`, `expressed_blood`, drop `is_purecontrol`) → subset to the 48 paired
   samples → gene collapse. The expression source pivoted to **GSE109597** (Affymetrix,
   84 samples, already RMA-normalized), so the Illumina flags/linkage don't apply — the
   rewrite parses the GSE series matrix and maps Affy probesets→genes (GPL570). See the
   pivot note + [[gse109597-expression-pivot]].
3. **`03_concordance.R`** — three-signal framework: (1) diet→methylation EWAS,
   (2) diet→expression DE, (3) methylation↔expression eQTM. Pair is *concordant* when
   `sign(b_diet_expr) == sign(b_diet_meth × b_eqtm)`. limma moderated stats; covariates
   Age+Sex+BMI, cell types adjust the methylation model (one dropped as reference to
   avoid collinearity — proportions sum to 1).
   **⚠️ Signal (3), the eQTM, assumed methylation + expression in the SAME samples —
   impossible now that expression is a separate cohort (GSE109597).** The framework
   must be redesigned (e.g. drop the within-subject eQTM; test cross-cohort gene-level
   direction agreement between diet→methylation and obesity→expression) — next-session
   work. `03a`/`03b`/`03c` (EWAS, hit annotation, candidate enrichment) are unaffected —
   they are methylation-only and already run.
4. **`04_dmr.R`** — DMRcate region-level analysis, companion to `03`'s per-CpG
   EWAS. Same `DietScore + Age + Sex + BMI + 5 cell-type` design → diet-associated
   DMRs + overlapping genes, Klemp-comparable params (λ=1000, C=2, >2 CpGs,
   min-smoothed-FDR<5%, |meandiff|≥2%). Scope: DMR calling + gene annotation only
   (no expression concordance). Flags DMRs overlapping the candidate list.
5. **Candidate-gene follow-up** — `data/raw/annotation/candidate_genes.txt`
   (25 genes), tagged in `03` for a powered arm. See "Candidate gene sourcing".
6. **Pathway enrichment** — clusterProfiler on concordant genes (not yet scripted).

Note: the genome-wide EWAS lives inside `03_concordance.R` (its stage 1) and is ALSO
available standalone as **`03a_ewas.R`** — Sections 0–3 carved out, sharing `03`'s
config + stats verbatim, so the methylation discovery scan can run before the
expression data arrives. The two produce the identical EWAS table (`03a` was added
this session because expression is still blocked).
Caveat carried in `04`: Klemp's 2% mean-diff floor was for a BINARY contrast; with
continuous `DietScore` `meandiff` is per-unit, so that floor is scale-dependent and
may over-filter (the script logs a warning if it removes all FDR-significant DMRs).

## Candidate gene sourcing & provenance
- **Sourced from the real paper, not memory.** The 25 genes were extracted directly
  from the Klemp et al. 2022 PDF (`data/raw/`) — Table 2 (top DMPs, p.6) and Table 3
  (top DMRs, p.10), read via rendered page images because the tables are rotated.
  Every CpG ID was cross-checked against machine-extracted text. Tooling: `pymupdf`
  (installed to the user's Python 3.13; poppler/pdftotext were unavailable).
- **Smoking-confound exclusion.** The 4 canonical smoking-methylation markers in
  those tables (`AHRR`, `F2RL3`, `GFI1`, `RARA`) were **deliberately excluded** — the
  paper flags them as its strongest effects (p<1e-11) and they are smoking- not
  diet-driven; the models already adjust for `Smoking_score`. This is the
  "metabolic-focused subset" (user's choice).
- **Caveat carried in-file:** these are Klemp's lifestyle-**composite** (healthy vs
  unhealthy) hits, used as an a-priori prior — NOT a diet-specific gene set. The
  candidate arm is hypothesis-generating.
- **Verified testable:** all 25 have ≥1 probe in `Gx_probe_info` (`symbol_INGENUITY`).
  Legacy alias `KIAA1026` (Table 2) was resolved to current symbol **`KAZN`** via the
  `synonyms_INGENUITY` field — using the legacy symbol would have silently dropped it.
- `03` parses the list stripping `#` comments (full-line and inline).

## Tech stack
- **R** for genomics statistics: `minfi`, `ChAMP`, `limma`, `DMRcate` (DMR calling),
  `sva`, `FlowSorted.Blood.EPIC`, `clusterProfiler`, `org.Hs.eg.db`, plus `maxprobes`
  (cross-reactive list, from GitHub — needs `minfiData`). Setup:
  `scripts/install_r_packages.R` — pins **Bioconductor 3.23**, installs
  non-interactively into the user library (`R_LIBS_USER`), verifies the objects the
  pipeline needs, and records `results/qc/r_package_versions.txt`.
  `scripts/install_flowsorted.R` is a targeted re-installer for the cell-type
  reference. **Host env:** R 4.6.1; `Rscript` is NOT on PATH in Git Bash
  (`C:\Program Files\R\R-4.6.1\bin\Rscript.exe`); ~7.8 GB RAM and segfault-prone under
  rapid `Rscript -e` calls, so run installs/checks as script files, not one-liners
  (see [[houseman-lowmem-and-host-limits]]).
- **Python** for data wrangling / glue / the eventual dashboard (`requirements.txt`:
  pandas, numpy, openpyxl, matplotlib, seaborn, jupyter).

## Repo layout
`data/raw/` (originals, gitignored) · `data/processed/` (normalized outputs) ·
`scripts/` (numbered pipeline) · `results/` (tables, figures, QC) · `notebooks/`
(exploratory) · `docs/` (writeups).

## End goal
Portfolio project → eventual **interactive web dashboard** presenting the concordance
results.

## Conventions
- Both assays key on methylation `Sample_Name`.
- Processed objects → `data/processed/*.rds`; QC → `results/qc/`; result tables →
  `results/tables/`; figures → `results/figures/`.
- Scripts run non-interactively (`Rscript scripts/NN_*.R`), self-contained `cfg` list
  at top, fail-fast file checks, and each writes a `sessionInfo` record.
- **Large outputs stay out of git:** `data/processed/*.rds` are gitignored; full
  genome-wide result tables (e.g. the 806k-row EWAS, `ewas_diet_methylation.csv.gz`)
  are gitignored + regenerable — only trimmed subsets are tracked (e.g.
  `ewas_diet_significant.csv`, p<1e-3 hits). Big scatter plots are rasterized PNGs,
  not multi-MB vector PDFs (an 806k-point volcano PDF was ~31 MB → PNG ~30 KB).
