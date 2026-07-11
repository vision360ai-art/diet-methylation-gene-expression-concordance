#!/usr/bin/env Rscript
# =============================================================================
# 01b_finalize_from_grset.R  --  recover the tail end of 01 from a saved grSet.
#
# 01's final normalization succeeded and grSet_normalized_filtered.rds was written,
# but the subsequent saveRDS(m_values) failed on a transient antivirus file lock
# ("error writing to connection"), leaving m_values.rds truncated and the cell-
# corrected pheno unwritten. Rather than re-run the crash-prone IDAT read +
# normalization, this reconstructs the remaining outputs deterministically from the
# saved grSet -- identical math to 01 sections 8-9, with a retrying writer.
#
# Idempotent: safe to re-run. Regenerates m_values.rds + pheno (rds/csv); also
# re-derives beta the same way 01 does and rewrites it so the whole set is
# guaranteed mutually consistent.
# =============================================================================

user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib) && dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages(library(minfi))

cfg <- list(out_dir = "data/processed")
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))
cellTypes <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")

grset_path <- file.path(cfg$out_dir, "grSet_normalized_filtered.rds")
if (!file.exists(grset_path)) stop("Missing ", grset_path, " -- re-run 01.")

log_msg("Loading %s ...", grset_path)
grSet <- readRDS(grset_path)
log_msg("grSet: %d probes x %d samples", nrow(grSet), ncol(grSet))

pheno <- as.data.frame(pData(grSet))
have_ct <- cellTypes %in% colnames(pheno)
if (!all(have_ct))
  stop("grSet pData is MISSING cell-type columns: ",
       paste(cellTypes[!have_ct], collapse = ", "),
       " -- the cell counts did not propagate; re-run 01 with the cache present.")
log_msg("Confirmed cell-type columns present. Means: %s",
        paste(sprintf("%s=%.3f", cellTypes, colMeans(pheno[, cellTypes])), collapse = ", "))

# ---- sections 8-9 of 01, verbatim math ------------------------------------
beta <- getBeta(grSet)
beta[beta < 1e-6]     <- 1e-6
beta[beta > 1 - 1e-6] <- 1 - 1e-6
mval <- log2(beta / (1 - beta))
log_msg("Beta matrix: %d x %d | M matrix: %d x %d",
        nrow(beta), ncol(beta), nrow(mval), ncol(mval))

robust_saveRDS <- function(obj, path, tries = 4, wait = 5) {
  for (i in seq_len(tries)) {
    ok <- tryCatch({ saveRDS(obj, path); TRUE },
                   error = function(e) { log_msg("  saveRDS(%s) attempt %d/%d failed: %s",
                                                 basename(path), i, tries, conditionMessage(e)); FALSE })
    if (ok) return(invisible(TRUE))
    if (file.exists(path)) unlink(path)
    Sys.sleep(wait); invisible(gc())
  }
  stop("Could not write ", path, " after ", tries, " attempts.")
}

robust_saveRDS(beta,  file.path(cfg$out_dir, "beta_values.rds"))
robust_saveRDS(mval,  file.path(cfg$out_dir, "m_values.rds"))
robust_saveRDS(pheno, file.path(cfg$out_dir, "pheno_methylation.rds"))
write.csv(pheno, file.path(cfg$out_dir, "pheno_methylation.csv"), row.names = FALSE)
log_msg("Wrote beta_values.rds, m_values.rds, pheno_methylation.{rds,csv}")
log_msg("DONE.")
