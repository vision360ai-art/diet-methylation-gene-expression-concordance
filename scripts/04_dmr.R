# =============================================================================
# 04_dmr.R
#
# Differentially methylated REGION (DMR) calling for diet quality, LIFE-Adult
# (n=48), using DMRcate. Region-level companion to the per-CpG EWAS in 03.
#
# WHY REGIONS
#   DMRs aggregate signal across neighbouring correlated CpGs, so a region can
#   reach significance where no single probe does -- useful at n=48. This is also
#   the method Klemp et al. 2022 used for their Table 3, so results are directly
#   comparable to the paper that seeded the candidate-gene list.
#
# SCOPE: DMR calling + gene annotation only. Expression concordance stays in 03
#   (per-CpG). This script produces a ranked table of diet-associated DMRs and the
#   genes they overlap; it does NOT test expression here.
#
# MODEL: same design as the EWAS in 03 --
#   M-value ~ DietScore + Age + Sex + BMI + 5 cell-type proportions (Neu dropped).
#
# PARAMETERS: matched to Klemp et al. for comparability (lambda=1000, C=2,
#   >2 CpGs, min-smoothed-FDR <5%, mean diff >=|2%|). See the DESIGN NOTE on the
#   mean-difference floor -- it is contrast-scale-dependent.
#
# Inputs  (from script 01): data/processed/m_values.rds, pheno_methylation.rds
# Outputs : results/tables/dmrs_diet_methylation.csv, data/processed/dmrs.rds,
#           results/figures/dmr_effect_vs_fdr.pdf
#
# STATUS: SKETCH. Cannot run until 01 has produced its processed objects (blocked
#   on the not-yet-downloaded IDATs).
#
# Run non-interactively:  Rscript scripts/04_dmr.R
# =============================================================================

suppressPackageStartupMessages({
  library(DMRcate)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(data.table)
})

# ----------------------------------------------------------------------------
# 0. Configuration
# ----------------------------------------------------------------------------
cfg <- list(
  proc_dir   = "data/processed",
  tbl_dir    = "results/tables",
  fig_dir    = "results/figures",

  # analysis variable + adjustment set (must match 03 for a coherent story)
  exposure   = "DietScore",
  covariates = c("Age", "Sex", "BMI"),
  cell_types    = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
  cell_ref_drop = "Neu",

  # DMRcate parameters (Klemp-comparable)
  arraytype   = "EPIC",   # DMRcate >=2.14 may require "EPICv1" for the 850K v1 array
  lambda      = 1000,     # Gaussian kernel bandwidth (bp)
  C           = 2,        # kernel scaling factor
  cpg_fdr     = 0.05,     # per-CpG FDR passed to cpg.annotate
  min_cpgs    = 2,        # keep regions with no.cpgs > this (Klemp: >2)
  dmr_fdr     = 0.05,     # region-level min-smoothed-FDR cutoff
  # DESIGN NOTE -- mean-difference floor is CONTRAST-SCALE-DEPENDENT.
  #   Klemp used |mean diff| >= 2% for a BINARY healthy-vs-unhealthy contrast.
  #   Here the exposure is CONTINUOUS DietScore, so `meandiff` is the change PER
  #   UNIT of DietScore (which spans ~21 units) -- a 0.02 per-unit floor is far
  #   stricter than Klemp's and will likely remove everything. Set to 0 to disable,
  #   or reason about it as effect across the DietScore range. Kept at 0.02 to
  #   honour "match Klemp", but the run logs how many DMRs it removes so the effect
  #   is visible, not silent.
  min_meandiff = 0.02,

  seed = 1234
)

set.seed(cfg$seed)
dir.create(cfg$tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$fig_dir, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))
rds <- function(f) { p <- file.path(cfg$proc_dir, f)
  if (!file.exists(p)) stop("Missing input: ", p, " -- run script 01 first.")
  readRDS(p) }

# small helper: pick the first present column name (DMRcate names vary by version)
pick_col <- function(df, candidates, what) {
  hit <- intersect(candidates, names(df))
  if (!length(hit)) stop(sprintf("No %s column found; looked for: %s. Got: %s",
                                 what, paste(candidates, collapse = ", "),
                                 paste(names(df), collapse = ", ")))
  hit[1]
}

# ----------------------------------------------------------------------------
# 1. Load methylation + build the (EWAS-matching) design
# ----------------------------------------------------------------------------
mval <- rds("m_values.rds")             # CpG x sample (M-values)
ph   <- rds("pheno_methylation.rds")    # DietScore, Age, Sex, BMI, cell props
ph   <- ph[colnames(mval), , drop = FALSE]
ph$Sex <- droplevels(factor(ph$Sex))

cells_used <- setdiff(cfg$cell_types, cfg$cell_ref_drop)
terms <- c(cfg$exposure, cfg$covariates, cells_used)
terms <- terms[terms %in% names(ph)]
design <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), data = ph)
mval   <- mval[, rownames(design), drop = FALSE]   # align samples to design rows
diet_col <- match(make.names(cfg$exposure), colnames(design))
log_msg("Design: ~ %s  (%d samples, DietScore = column %d)",
        paste(colnames(design)[-1], collapse = " + "), nrow(design), diet_col)

# ----------------------------------------------------------------------------
# 2. Per-CpG annotation (limma inside DMRcate) on the DietScore coefficient
# ----------------------------------------------------------------------------
log_msg("cpg.annotate over %d CpGs (arraytype = %s)...", nrow(mval), cfg$arraytype)
annot <- cpg.annotate(
  object        = mval,
  datatype      = "array",
  what          = "M",
  arraytype     = cfg$arraytype,
  analysis.type = "differential",
  design        = design,
  coef          = diet_col,
  fdr           = cfg$cpg_fdr
)

# ----------------------------------------------------------------------------
# 3. Agglomerate CpGs into regions + annotate genes (hg19)
# ----------------------------------------------------------------------------
log_msg("dmrcate (lambda = %d, C = %d)...", cfg$lambda, cfg$C)
dmrs <- dmrcate(annot, lambda = cfg$lambda, C = cfg$C)
ranges <- extractRanges(dmrs, genome = "hg19")
saveRDS(ranges, file.path(cfg$proc_dir, "dmrs.rds"))

df <- as.data.frame(ranges)
fdr_col  <- pick_col(df, c("min_smoothed_fdr", "HMFDR", "Stouffer", "Fisher"), "FDR")
diff_col <- pick_col(df, c("meandiff", "meanbetafc"), "mean-difference")
ncpg_col <- pick_col(df, c("no.cpgs", "no.probes"), "CpG-count")
log_msg("Raw DMRs: %d (FDR col = '%s', diff col = '%s').",
        nrow(df), fdr_col, diff_col)

# ----------------------------------------------------------------------------
# 4. Filter to significant regions (Klemp-comparable)
# ----------------------------------------------------------------------------
n0 <- nrow(df)
keep_cpgs <- df[[ncpg_col]] >  cfg$min_cpgs
keep_fdr  <- df[[fdr_col]]  <  cfg$dmr_fdr
keep_diff <- abs(df[[diff_col]]) >= cfg$min_meandiff
keep_all  <- keep_cpgs & keep_fdr & keep_diff
log_msg("Filter counts: >%d CpGs %d | FDR<%.2f %d | |meandiff|>=%.3f %d (of %d)",
        cfg$min_cpgs, sum(keep_cpgs), cfg$dmr_fdr, sum(keep_fdr),
        cfg$min_meandiff, sum(keep_diff), n0)
if (cfg$min_meandiff > 0 && sum(keep_cpgs & keep_fdr) > 0 && sum(keep_all) == 0) {
  log_msg("WARNING: the mean-diff floor removed ALL FDR-significant DMRs -- likely the")
  log_msg("         continuous-exposure scale issue (see cfg DESIGN NOTE). Try min_meandiff=0.")
}

sig <- df[keep_all, , drop = FALSE]
sig$direction <- ifelse(sig[[diff_col]] < 0, "hypomethylated", "hypermethylated")
sig <- sig[order(sig[[fdr_col]]), ]
log_msg("Significant DMRs: %d (%d hypo, %d hyper).",
        nrow(sig), sum(sig$direction == "hypomethylated"),
        sum(sig$direction == "hypermethylated"))

# ----------------------------------------------------------------------------
# 5. Gene annotation + candidate-list cross-reference (annotation only)
# ----------------------------------------------------------------------------
gene_col <- pick_col(sig, c("overlapping.genes", "overlapping.promoters"), "gene")
cand_file <- "data/raw/annotation/candidate_genes.txt"
if (nrow(sig) && file.exists(cand_file)) {
  cand <- sub("#.*$", "", readLines(cand_file)); cand <- unique(trimws(cand))
  cand <- cand[cand != ""]
  sig$in_candidate_list <- vapply(strsplit(sig[[gene_col]], ",\\s*"),
    function(gs) any(trimws(gs) %in% cand), logical(1))
  log_msg("DMRs overlapping a candidate gene: %d.", sum(sig$in_candidate_list, na.rm = TRUE))
}

fwrite(sig, file.path(cfg$tbl_dir, "dmrs_diet_methylation.csv"))

# ----------------------------------------------------------------------------
# 6. Figure + session record
# ----------------------------------------------------------------------------
if (nrow(df)) {
  pdf(file.path(cfg$fig_dir, "dmr_effect_vs_fdr.pdf"), width = 6, height = 6)
  plot(df[[diff_col]], -log10(pmax(df[[fdr_col]], 1e-300)), pch = 16, cex = 0.5,
       col = ifelse(keep_all, "firebrick", "grey60"),
       xlab = "Region mean methylation difference", ylab = "-log10 min-smoothed FDR",
       main = "Diet-associated DMRs (red = retained)")
  abline(v = c(-cfg$min_meandiff, cfg$min_meandiff), h = -log10(cfg$dmr_fdr),
         col = "grey80", lty = 2)
  dev.off()
}
writeLines(capture.output(sessionInfo()),
           file.path(cfg$tbl_dir, "sessionInfo_04_dmr.txt"))
log_msg("Done. DMR table in %s/, GRanges in %s/dmrs.rds.", cfg$tbl_dir, cfg$proc_dir)
