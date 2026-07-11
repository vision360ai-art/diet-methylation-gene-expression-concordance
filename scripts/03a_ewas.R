#!/usr/bin/env Rscript
# =============================================================================
# 03a_ewas.R  --  standalone diet -> methylation EWAS.
#
# This is Sections 0-3 of 03_concordance.R carved out so the methylation-side
# discovery scan can run NOW, independent of the expression data (the values
# matrix / expr_gene_paired.rds is not yet available). It intentionally shares
# 03's config values and statistical recipe VERBATIM, so when 03 later runs its
# Section 3 it reproduces this table bit-for-bit.
#
# Scope: (1) diet -> methylation EWAS only. It does NOT do CpG->gene mapping,
#        diet->expression DE, eQTM, or concordance -- those need expression data
#        and live in 03_concordance.R.
#
# Model:  M-value ~ DietScore + Age + Sex + BMI + 5 cell-type proportions
#         (6 Houseman proportions from 01; Neu dropped as reference -- they sum to
#          ~1, so keeping all 6 is collinear with the intercept). limma moderated t.
#
# Inputs  : data/processed/m_values.rds        (CpG x sample, M-values)
#           data/processed/pheno_methylation.rds (DietScore, Age, Sex, BMI + 6 cells)
# Outputs : results/tables/ewas_diet_significant.csv  (TRACKED: nominal p<1e-3 hits,
#           ~1k rows -- the discovery pool that feeds 03's concordance stage)
#           results/tables/ewas_diet_methylation.csv.gz  (GITIGNORED: full 806k-row
#           table, ~35 MB gzipped; regenerate via this script. fread() auto-gunzips)
#           results/figures/ewas_volcano.png   (rasterized: 806k vector points -> 31 MB
#           PDF; PNG is a few hundred KB and overplots identically)
#           results/tables/sessionInfo_03a_ewas.txt
#
# Run:  Rscript scripts/03a_ewas.R
# =============================================================================

user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib) && dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(limma)
  library(data.table)
})

# ---- config (kept in sync with 03_concordance.R) --------------------------
cfg <- list(
  proc_dir   = "data/processed",
  tbl_dir    = "results/tables",
  fig_dir    = "results/figures",
  exposure   = "DietScore",
  covariates = c("Age", "Sex", "BMI"),
  cell_types      = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
  cell_ref_drop   = "Neu",
  ewas_p_discovery = 1e-3,   # only used for the volcano threshold line here
  seed = 1234
)
set.seed(cfg$seed)
dir.create(cfg$tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$fig_dir, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))
rds <- function(f) { p <- file.path(cfg$proc_dir, f)
  if (!file.exists(p)) stop("Missing input: ", p, " -- run 01 (and 01a) first.")
  readRDS(p) }

# ---- 1. load & harmonize (methylation only; no expression join) ----------
mval       <- rds("m_values.rds")
pheno_meth <- rds("pheno_methylation.rds")

# EWAS-only "common" set = methylation samples that also have a phenotype row.
common <- intersect(colnames(mval), rownames(pheno_meth))
log_msg("Methylation samples: %d; with phenotype: %d -> EWAS on %d.",
        ncol(mval), nrow(pheno_meth), length(common))
if (length(common) < 10) stop("Too few samples (", length(common), ").")
mval <- mval[, common, drop = FALSE]
ph   <- pheno_meth[common, , drop = FALSE]
ph$Sex <- droplevels(factor(ph$Sex))

# fail loudly if the cell-type covariates aren't present (i.e. Houseman skipped)
missing_cells <- setdiff(cfg$cell_types, names(ph))
if (length(missing_cells))
  stop("pheno lacks cell-type columns: ", paste(missing_cells, collapse = ", "),
       " -- run 01a_cell_counts.R + 01 so the correction is in pheno.")

# ---- 2. methylation design matrix (identical to 03's build_design) --------
cells_used <- setdiff(cfg$cell_types, cfg$cell_ref_drop)
terms <- c(cfg$exposure, cfg$covariates, cells_used)
terms <- terms[terms %in% names(ph)]
mm    <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), data = ph)
coef  <- make.names(cfg$exposure)
log_msg("Methylation model:  ~ %s", paste(colnames(mm)[-1], collapse = " + "))
if (nrow(mm) < length(common))
  log_msg("WARNING: %d samples dropped by model.matrix (missing covariates).",
          length(common) - nrow(mm))

# ---- 3. diet -> methylation EWAS (limma moderated t) ---------------------
log_msg("EWAS: M-value ~ diet + covars + cell-type over %d CpGs, n=%d...",
        nrow(mval), nrow(mm))
fit_m  <- eBayes(lmFit(mval[, rownames(mm)], mm))
ewas   <- topTable(fit_m, coef = coef, number = Inf, sort.by = "p")
ewas$CpG <- rownames(ewas)
# Full table: gzip on write (fwrite auto-detects .gz) -- ~35 MB, GITIGNORED (too big
# to track; regenerable). Read back with data.table::fread(), which auto-decompresses.
fwrite(ewas, file.path(cfg$tbl_dir, "ewas_diet_methylation.csv.gz"))
# Tracked artifact: just the nominal-significant discovery pool (small).
sig <- ewas[ewas$P.Value < cfg$ewas_p_discovery, ]
sig <- sig[order(sig$P.Value), ]
fwrite(sig, file.path(cfg$tbl_dir, "ewas_diet_significant.csv"))
log_msg("Wrote %d significant CpGs (p<%.0e) -> ewas_diet_significant.csv (tracked).",
        nrow(sig), cfg$ewas_p_discovery)
log_msg("EWAS done. Genome-wide FDR<0.05: %d CpGs (expected ~0 at this n).",
        sum(ewas$adj.P.Val < 0.05))
log_msg("Nominal p<%.0e: %d CpGs; p<1e-4: %d; p<1e-5: %d.",
        cfg$ewas_p_discovery, sum(ewas$P.Value < cfg$ewas_p_discovery),
        sum(ewas$P.Value < 1e-4), sum(ewas$P.Value < 1e-5))

# top-hits preview to the log
top <- head(ewas[order(ewas$P.Value), c("CpG","logFC","P.Value","adj.P.Val")], 8)
cat("Top CpGs (diet->methylation):\n"); print(top, row.names = FALSE)

# ---- volcano (rasterized PNG: 806k points overplot the same, ~KB not ~31 MB) ----
png(file.path(cfg$fig_dir, "ewas_volcano.png"), width = 1200, height = 1200, res = 200)
plot(ewas$logFC, -log10(ewas$P.Value), pch = 16, cex = 0.3, col = "grey60",
     xlab = "Diet effect on M-value", ylab = "-log10 p",
     main = "EWAS: diet -> methylation")
abline(h = -log10(cfg$ewas_p_discovery), col = "red", lty = 2)
dev.off()

writeLines(capture.output(sessionInfo()),
           file.path(cfg$tbl_dir, "sessionInfo_03a_ewas.txt"))
log_msg("Done. Table: %s/ewas_diet_methylation.csv.gz | Figure: %s/ewas_volcano.png",
        cfg$tbl_dir, cfg$fig_dir)
