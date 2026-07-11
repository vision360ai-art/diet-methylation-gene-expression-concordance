#!/usr/bin/env Rscript
# =============================================================================
# 01a_cell_counts.R  --  Houseman blood cell-type deconvolution, LOW-MEMORY.
#
# Why this is a separate, hand-rolled script instead of a plain estimateCellCounts2
# call inside 01:
#   estimateCellCounts2() co-normalizes the user samples together with the full
#   ~37-sample sorted-blood reference (a combined ~73 x 1.05M-probe object,
#   preprocessNoob). Its peak memory is ~5.5-6 GB. This host has 7.8 GB RAM but
#   only ~5.7 GB of commit heaadroom for a single process once the OS + apps are
#   loaded, so the combined object OOMs -- and batching the USER samples does not
#   help, because the ~37 reference samples dominate the combined object.
#
#   This script uses the SAME method via FlowSorted.Blood.EPIC's exported low-level
#   API, which never materialises the combined object:
#     1. preprocessNoob() on the USER RGChannelSet ALONE  (~2 GB peak, fits).
#        noob is a WITHIN-array background correction (each array uses its own
#        out-of-band probes; no cross-sample borrowing), so user-only noob yields
#        beta values IDENTICAL to what co-normalisation would produce for them.
#     2. IDOLOptimizedCpGs.compTable -- the PRECOMPUTED reference cell-type mean
#        methylation at the IDOL library CpGs (this is exactly the coefficient
#        matrix estimateCellCounts2 builds internally from the reference).
#     3. projectCellType_CP() -- the constrained-projection (Houseman) solver,
#        run per sample against that reference. This is the same function
#        estimateCellCounts2 calls for the final projection.
#   Result is statistically equivalent to estimateCellCounts2(processMethod=
#   "preprocessNoob", probeSelect="IDOL") at a fraction of the peak memory.
#
# It also caches the result: 01 reads data/processed/cell_counts.csv if present
# instead of recomputing, so 01 re-runs are fast and never re-hit the memory wall.
#
# Input  : data/processed/rgSet_raw.rds   (written by 01, sampleNames = Sample_Name)
# Output : data/processed/cell_counts.csv (rows = Sample_Name, 6 cell proportions)
#          results/qc/qc_cell_composition.pdf
#
# Run:  Rscript scripts/01a_cell_counts.R   (AFTER 01 has written rgSet_raw.rds)
# =============================================================================

# Target the user-writable library first so FlowSorted.Blood.EPIC resolves.
user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib) && dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(minfi)
  library(FlowSorted.Blood.EPIC)
})

cfg <- list(
  rgset     = "data/processed/rgSet_raw.rds",
  out       = "data/processed/cell_counts.csv",
  qc_pdf    = "results/qc/qc_cell_composition.pdf",
  qc_dir    = "results/qc",
  cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
)
dir.create(cfg$qc_dir, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))

if (!file.exists(cfg$rgset))
  stop("Missing ", cfg$rgset, " -- run 01_preprocess_methylation.R first.")

# Serial read/preprocess: minfi's BiocParallel default (SnowParam) spawns workers
# that each load the ~300 MB EPIC manifest -> OOM. Same rationale as 01.
BiocParallel::register(BiocParallel::SerialParam())

# ---- 1. user-only noob ----------------------------------------------------
log_msg("Loading raw RGChannelSet: %s", cfg$rgset)
rgSet <- readRDS(cfg$rgset)
log_msg("RGChannelSet: %d probes x %d samples", nrow(rgSet), ncol(rgSet))

log_msg("preprocessNoob() on user samples alone...")
mset <- preprocessNoob(rgSet)
rm(rgSet); invisible(gc())
beta <- getBeta(mset)
rm(mset); invisible(gc())
log_msg("User beta matrix: %d probes x %d samples", nrow(beta), ncol(beta))

# ---- 2. precomputed IDOL reference coefficient matrix ---------------------
# IDOLOptimizedCpGs.compTable rows = IDOL library CpGs; it carries the reference
# cell-type mean methylation columns (the projection coefficients).
comp <- get("IDOLOptimizedCpGs.compTable")
if (!all(cfg$cellTypes %in% colnames(comp)))
  stop("compTable missing expected cell-type columns. Present: ",
       paste(colnames(comp), collapse = ", "))
coefs <- as.matrix(comp[, cfg$cellTypes, drop = FALSE])

probes <- intersect(rownames(coefs), rownames(beta))
log_msg("IDOL library CpGs: %d in compTable, %d present in user data (using %d).",
        nrow(coefs), length(probes), length(probes))
if (length(probes) < 0.8 * nrow(coefs))
  log_msg("  WARNING: >20%% of IDOL CpGs absent from user data -- proportions may be degraded.")

# ---- 3. constrained-projection (Houseman) per sample ----------------------
log_msg("Projecting cell-type proportions (projectCellType_CP)...")
counts <- projectCellType_CP(
  beta[probes, , drop = FALSE],
  coefs[probes, , drop = FALSE],
  nonnegative = TRUE,
  lessThanOne = FALSE
)
cell_counts <- as.data.frame(counts)
cell_counts <- cell_counts[colnames(beta), cfg$cellTypes, drop = FALSE]  # order

# ---- 4. write cache + QC --------------------------------------------------
write.csv(cell_counts, cfg$out, row.names = TRUE)
log_msg("Wrote %d x %d cell-proportion table -> %s",
        nrow(cell_counts), ncol(cell_counts), cfg$out)
log_msg("Cell proportions (mean): %s",
        paste(sprintf("%s=%.3f", colnames(cell_counts), colMeans(cell_counts)),
              collapse = ", "))

pdf(cfg$qc_pdf, width = 10, height = 5)
boxplot(cell_counts, ylab = "Estimated proportion", las = 2,
        main = "Houseman blood cell-type proportions (noob + IDOL projection)",
        col = "grey80")
dev.off()
log_msg("Wrote QC plot -> %s", cfg$qc_pdf)
