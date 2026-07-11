# Project Concept & Execution Plan

Plain-language, start-to-finish plan for the diet–methylation–expression
concordance project. Doubles as the portfolio narrative / interview talking
points. For the technical brief (data quirks, models, conventions) see
[`../CLAUDE.md`](../CLAUDE.md).

**Status legend:** ✅ done · 🚧 in progress · ⏳ blocked / waiting · ⬜ not started

---

## Where things stand (2026-07-10)

**Phase: pipeline fully scaffolded, blocked on data.** Every analysis script is
written and committed; nothing has been *run* yet because two raw-data files
haven't arrived.

- ✅ **Code complete (pending real data):** preprocessing + analysis scripts
  `01`→`03` written, committed, and pushed; R environment setup
  (`install_r_packages.R`); technical brief (`CLAUDE.md`) and this plan.
- ✅ **Candidate gene list built from source:** 25 genes extracted from the real
  Klemp et al. PDF, verified against the expression annotation, committed.
- ⏳ **Blocking everything downstream — two files not yet in repo:**
  methylation **IDATs** (0 / 98 present) and the **expression values matrix**
  (missing). The sample sheet, sample info, and probe annotation *are* in hand.
- ⬜ **No analysis has run:** `data/processed/` and `results/` are empty; there are
  **no statistical results yet** — expected at this stage, not a gap.
- ▶︎ **Immediate next action:** obtain the IDATs + expression matrix, then execute
  `01` → `02` → `03` and review QC before trusting any output.

---

## 1. The question
Using real data from 48 people (blood samples, LIFE-Adult cohort, Germany),
determine whether diet quality relates to genes getting chemically "muted"
(methylation) in a way that actually shows up as the gene being less active
(expression) — not just methylation changes that happen but don't do anything.

## 2. Get the raw data ready 🚧 *(the current bottleneck)*
- ✅ Sample sheets and phenotype info — done, verified.
- ✅ Candidate gene list from prior research on this cohort (Klemp et al.) —
  **extracted directly from the real PDF** (Tables 2 & 3), 25 genes, verified
  against the expression annotation. See step 7.
- ⏳ Raw methylation chip files (IDATs) — not yet in repo; expected in
  `data/raw/idat/`.
- ⏳ Raw gene expression values matrix — not yet in repo. (The probe annotation
  and sample info are already in hand.)

## 3. Clean and prepare the methylation data ✅ *written* · ⏳ *not yet run (needs IDATs)*
`scripts/01_preprocess_methylation.R`
- Load raw chip signal.
- Quality-check each sample, drop anything broken (detection p-values, sex check).
- Estimate each person's blood cell-type mixture (Houseman), so diet's effect on
  immune-cell composition isn't mistaken for a direct effect on gene methylation.
- Normalize to remove machine/technical noise while preserving real biology
  (functional normalization).
- Filter out unreliable measurement spots.
- Save a clean, ready-to-analyze dataset (Beta/M matrices + phenotype).

## 4. Clean and prepare the expression data ✅ *drafted* · ⏳ *waiting on values file*
`scripts/02_preprocess_expression.R`
- Load the already-normalized expression values (the LIFE-Adult team pre-processed
  this centrally).
- Apply their built-in quality flags.
- Map each measurement to the actual gene it represents.
- Subset down to just the 48 people who also have methylation data
  (linkage: expression `sampleID` == methylation `Synonym`).
- Save a clean, ready-to-analyze dataset.
- *Confirmed:* `probe_annot` wired to the real `Gx_probe_info` file. Still needs
  the per-sample values matrix to finalize.

## 5. Main statistical test + concordance ✅ *written* · ⏳ *needs real data to run*
`scripts/03_concordance.R` — **the core deliverable.**

> **Design decision:** the genome-wide test and the concordance analysis live in
> **one script (`03`)**, not two. (An earlier draft of this plan split them; we
> keep them together.) `03` runs a **three-signal** framework:

1. **Diet → methylation (genome-wide EWAS).** For every methylation site, does it
   move with `DietScore` after accounting for age, sex, BMI, and blood cell mix?
   Multiple-testing (FDR) correction is applied. Results are reported honestly
   even if little or nothing survives strict correction at n=48 — that is a
   legitimate, expected outcome at this sample size, not a failure.
2. **Diet → expression.** Does the corresponding gene's expression move with diet?
3. **Methylation ↔ expression (eQTM).** Is there a local link between the site's
   methylation and the gene's expression?

A CpG–gene pair is **concordant** when the diet→expression effect matches what the
diet→methylation effect predicts through the local eQTM link
(`sign(b_diet_expr) == sign(b_diet_meth × b_eqtm)`). Output: a ranked table of
functionally-backed diet-methylation links (concordant) vs. methylation changes
with no downstream effect (discordant).

## 6. Candidate-gene follow-up ✅ *gene list done* · ⏳ *runs with step 5*
`data/raw/annotation/candidate_genes.txt` (25 genes)
- Re-run a smaller, targeted version of the same test on genes already implicated
  in this cohort's prior published research (Klemp et al.), scoped to a
  **metabolic-focused subset** (the four canonical smoking-methylation markers
  were excluded, since the models already adjust for smoking).
- Much smaller multiple-testing penalty → better chance of a real, statistically
  solid hit even with only 48 people.
- *Caveat carried in the file:* these are lifestyle-**composite** hits used as an
  a-priori prior, not a diet-specific gene set — the candidate arm is
  hypothesis-generating.

## 7. Interpretation ⬜
- Pathway enrichment (`clusterProfiler`) on whatever gene list comes out of steps
  5–6 — turn a list of gene names into a biological story ("these cluster around
  inflammation," etc.).
- Write up findings honestly, including the sample-size limitation, as part of the
  portfolio narrative.

## 8. Portfolio polish 🚧
- ✅ Groundwork: GitHub repo live, MIT-licensed, documented; `CLAUDE.md` in place
  for session continuity; git identity configured.
- ⬜ Clean README, key plots, plain-language summary of what was found (doubles as
  interview talking points).

## 9. Web dashboard ⬜ *(goal stated, no build decisions yet)*
Once real results exist, build an interactive, polished dashboard presenting the
concordance findings. Tech stack (Streamlit / Dash / etc.) to be decided once we
know the actual shape of the output data.

---

### DMR-level arm ✅ *written* · ⏳ *needs real data to run*
`scripts/04_dmr.R` — region-level companion to the per-CpG EWAS, using `DMRcate`
with Klemp-comparable parameters. Produces a ranked table of diet-associated DMRs
and their overlapping genes (scope: calling + annotation only; no expression
concordance). Note: Klemp's 2% mean-difference floor was calibrated for a binary
healthy-vs-unhealthy contrast; with continuous `DietScore` the region `meandiff` is
per-unit, so that floor is scale-dependent — the script warns if it over-filters.
