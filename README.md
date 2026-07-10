# Diet, DNA Methylation, and Gene Expression Concordance in LIFE-Adult

## Research Question
In the LIFE-Adult cohort (n=48, blood, paired methylation + expression arrays), which CpG
methylation changes associated with diet quality (continuous DietScore) show concordant
changes in expression of the corresponding gene — independent of BMI, age, and sex?

## Background
This project applies an established methylation-expression concordance framework (used
previously in cognition, grip strength, and disease-biomarker studies) to a dataset/question
combination that does not appear to have been previously published: the LIFE-Adult cohort's
diet-quality score, cross-referenced against its paired methylation and expression array data.

The original LIFE-Adult lifestyle-methylation paper (Klemp et al. 2022, *Clinical and
Translational Medicine*) analyzed methylation only. This project extends that work by
integrating the cohort's expression array data, which was collected on an overlapping
but not previously cross-referenced subset of subjects.

## Data Source
- **Cohort:** LIFE-Adult study, Leipzig, Germany (~10,000 participants overall)
- **Methylation:** Illumina Infinium MethylationEPIC BeadChip, blood, n=98 subjects
  (extreme healthy/unhealthy lifestyle groups)
- **Expression:** Illumina HumanHT-12 v4 Expression BeadChip, blood, n=1,883 subjects
- **Paired overlap:** n=48 subjects have both assays
- **Access:** Leipzig Health Atlas, https://www.health-atlas.de/studies/57 (open access)
- **Key variable:** `DietScore` (continuous), plus `PA_Score`, `Smoking_score`, `AlcScore`,
  and composite `HLS` (Healthy Lifestyle Score)

## Known Limitations (stated upfront, not discovered later)
- n=48 is underpowered for genome-wide FDR-corrected discovery across ~850K methylation
  probes — see Analysis Plan for how this is handled (genome-wide + candidate-gene approach)
- Methylation and expression are both array-based (not WGBS/RRBS or RNA-seq) — this is
  standard practice for this type of paired human cohort data, not a compromise
- Blood cell-type composition is a known confounder for blood methylation and is corrected
  for using the Houseman method (via `minfi`)

## Analysis Plan
1. **Preprocessing** — QC, normalization, batch correction (see `scripts/01_preprocess_methylation.R`,
   `scripts/02_preprocess_expression.R`)
2. **Cell-type correction** — Houseman deconvolution for methylation data
3. **Genome-wide association** — `limma` linear models, DietScore ~ methylation / expression,
   adjusted for age, sex, BMI, cell type; BH-FDR correction
4. **Concordance analysis** — map significant/top-ranked CpGs to genes, cross-reference against
   expression associations, classify concordant vs. discordant
5. **Candidate-gene follow-up** — targeted analysis on ~20-50 inflammation/metabolic genes
   (TNF, IL6, NFKB1, CRP, PPARG, LEP, ADIPOQ, FTO) to address multiple-testing power limits
6. **Interpretation** — KEGG/GO pathway enrichment on concordant gene set

## Repo Structure
```
data/raw/          - original downloaded files (not tracked in git, see .gitignore)
data/processed/    - cleaned/normalized data outputs
scripts/           - R and Python analysis scripts, numbered by pipeline order
results/           - figures, tables, output statistics
notebooks/         - exploratory analysis
docs/              - project documentation, writeups
```

## Requirements
See `scripts/install_r_packages.R` and `requirements.txt`

## Status
🚧 In progress — data sourcing and environment setup complete, preprocessing pipeline next.
