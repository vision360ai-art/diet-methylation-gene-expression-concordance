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
| `Gx_DataDictionary-21-06-15_v3.xlsx` | Column docs for sample_info + probe annotation | — | — |

- **Cross-assay linkage:** methylation `Synonym` column == expression `sampleID`.
  Verified overlap = **48 subjects** with both assays. Paired cohort is 25F / 23M
  (so `Sex` is a real covariate). `DietScore` range 3–24, continuous.
- Preprocessing relabels expression columns from `sampleID` → methylation
  `Sample_Name` so both assays share one sample key downstream.

## Known data quirks
- **Duplicate `HLS` column** in `SampleSheet_V2.csv`: appears twice — categorical
  (`low`/`high`) and numeric (composite Healthy Lifestyle Score). Load with
  `read.csv(check.names = FALSE)` and rename explicitly to **`HLS_group`** and
  **`HLS_score`** (done in `01_preprocess_methylation.R`).
- **IDAT files not yet in repo.** minfi needs the raw `*_Grn.idat`/`*_Red.idat`
  pairs; they are gitignored and expected in **`data/raw/idat/`**. Script 01 fails
  fast with a clear message if absent.
- **Expression values matrix not yet downloaded.** Data dictionary implies it ships
  already normalized/batch-corrected/imputed (with per-probe QC flags), so script 02
  loads a processed matrix + applies provider flags rather than doing bead-level
  normalization. Unknown filenames/columns are marked `# CONFIRM` in `02`'s config.

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
4. **Candidate-gene follow-up** — ~20–50 a-priori inflammation/metabolic genes
   (`data/raw/annotation/candidate_genes.txt`), tagged in `03` for a powered arm.
5. **Pathway enrichment** — clusterProfiler on concordant genes (not yet scripted).

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
