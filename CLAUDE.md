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

## Data structure (`data/raw/`, gitignored)
| File | What | n | Platform |
|---|---|---|---|
| `SampleSheet_V2.csv` | Methylation sample sheet + phenotypes | 98 | Illumina MethylationEPIC (~850K), blood |
| `Gx_sample_info-21-06-15_v3.txt` | Expression sample metadata | 1,883 | Illumina HumanHT-12 v4, blood |
| `Gx_probe_info-21-06-15_v3.txt` | Expression probe annotation + QC flags (48,106 probes × 42 cols) | — | — |
| `Gx_DataDictionary-21-06-15_v3.xlsx` | Column docs for sample_info + probe annotation | — | — |
| `Clinical…Klemp…2022.pdf` | Source paper (candidate-gene provenance, see below) | — | — |

- **Cross-assay linkage:** methylation `Synonym` column == expression `sampleID`.
  Verified overlap = **48 subjects** with both assays. Paired cohort is 25F / 23M
  (so `Sex` is a real covariate). `DietScore` range 3–24, continuous.
- Preprocessing relabels expression columns from `sampleID` → methylation
  `Sample_Name` so both assays share one sample key downstream.
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
- **Expression values matrix not yet downloaded.** The probe *annotation*
  (`Gx_probe_info`) and *sample* metadata ARE in hand, but the per-sample values
  matrix is not. Data dictionary implies it ships already normalized/batch-
  corrected/imputed (per-probe QC flags), so script 02 loads a processed matrix +
  applies provider flags rather than doing bead-level normalization. `02`'s
  `probe_annot` path is confirmed to `Gx_probe_info`; the `expr_matrix` filename/
  columns are still marked `# CONFIRM`.

## Current status & blockers (2026-07-10)
Pipeline is **fully scaffolded and pushed** (`01`→`03`, install script, candidate
list, docs) but **nothing has been run** — `data/processed/` and `results/` are
empty. **Two missing raw files block everything downstream:** the methylation
**IDATs** (0/98 present) and the **expression values matrix**. Full plain-language
status lives in [`docs/execution_plan.md`](docs/execution_plan.md); this file is the
technical brief.

## Known limitations (stated upfront, not compromises)
- **n=48 paired is underpowered** for genome-wide FDR across ~850K probes.
  Handled by treating EWAS as *discovery* (nominal-p + top-N) and controlling FDR
  *within* the CpG–gene pair set, plus a candidate-gene arm.
- **Both platforms are arrays** (EPIC + HT-12v4), not WGBS/RRBS or RNA-seq. This is
  standard practice for paired human-cohort data of this type, not a shortcut.
- **Blood cell composition** is the dominant methylation confounder — corrected via
  Houseman deconvolution (see pipeline).

## Analysis plan / pipeline (`scripts/`, numbered by order)
1. **`01_preprocess_methylation.R`** — minfi: read IDATs → QC (detP, sex check) →
   drop failed samples → **Houseman cell-type deconvolution** (`estimateCellCounts2`,
   FlowSorted.Blood.EPIC/IDOL; 6 proportions written into pheno) → functional
   normalization → probe filtering → Beta/M matrices.
2. **`02_preprocess_expression.R`** — load processed HT-12 matrix → provider QC-flag
   filtering (`probe_QCok`, `expressed_blood`, drop `is_purecontrol`) → subset to the
   48 paired samples → optional gene collapse. ComBat wired but OFF (data pre-corrected).
3. **`03_concordance.R`** — three-signal framework: (1) diet→methylation EWAS,
   (2) diet→expression DE, (3) methylation↔expression eQTM. Pair is *concordant* when
   `sign(b_diet_expr) == sign(b_diet_meth × b_eqtm)`. limma moderated stats; covariates
   Age+Sex+BMI, cell types adjust the methylation model (one dropped as reference to
   avoid collinearity — proportions sum to 1).
4. **`04_dmr.R`** — DMRcate region-level analysis, companion to `03`'s per-CpG
   EWAS. Same `DietScore + Age + Sex + BMI + 5 cell-type` design → diet-associated
   DMRs + overlapping genes, Klemp-comparable params (λ=1000, C=2, >2 CpGs,
   min-smoothed-FDR<5%, |meandiff|≥2%). Scope: DMR calling + gene annotation only
   (no expression concordance). Flags DMRs overlapping the candidate list.
5. **Candidate-gene follow-up** — `data/raw/annotation/candidate_genes.txt`
   (25 genes), tagged in `03` for a powered arm. See "Candidate gene sourcing".
6. **Pathway enrichment** — clusterProfiler on concordant genes (not yet scripted).

Note: the genome-wide EWAS lives **inside `03_concordance.R`** (its stage 1), not a
separate script — a deliberate decision (an earlier plan draft split them).
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
  `sva`, `FlowSorted.Blood.EPIC`, `clusterProfiler`, `org.Hs.eg.db`.
  Setup: `scripts/install_r_packages.R`.
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
