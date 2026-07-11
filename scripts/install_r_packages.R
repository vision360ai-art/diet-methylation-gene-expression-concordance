#!/usr/bin/env Rscript
# =============================================================================
# install_r_packages.R  --  reproducible R environment setup for the
# diet / methylation / expression concordance project.
#
# Reproducibility notes (why this is more than a bare install list):
#   * Installs into the USER-writable library (R_LIBS_USER), never Program Files
#     (which needs admin). We PREPEND it to .libPaths() so base/site packages stay
#     readable while every new install lands in the user library.
#   * PINS the Bioconductor release so the same package versions resolve on a fresh
#     machine (Bioc 3.23 <-> R 4.6.x). This is what makes "the final run" repeatable
#     rather than drifting with whatever Bioc is current.
#   * Installs non-interactively: update=FALSE + ask=FALSE (no "update all/some/none?"
#     prompt that stalls or errors under Rscript) and INSTALL_opts="--no-lock"
#     (staged-install locking otherwise fails on this Windows host -- this is the
#     flag whose absence made the earlier FlowSorted.Blood.EPIC install fail).
#   * VERIFIES at the end that every critical package loads AND that the specific
#     FlowSorted objects the memory-safe Houseman path (01a_cell_counts.R) uses
#     -- projectCellType_CP + IDOLOptimizedCpGs.compTable -- are present. A silent
#     partial install would otherwise only surface mid-pipeline.
#   * Writes a version manifest to results/qc/ so the exact environment is recorded.
#
# Pinned environment (captured 2026-07-10): R 4.6.1, Bioc 3.23,
#   FlowSorted.Blood.EPIC 2.16.0, minfi 1.58.0, limma 3.68.4, DMRcate 3.8.0,
#   sva 3.60.0, IlluminaHumanMethylationEPICmanifest 0.3.0.
#
# Run:  Rscript scripts/install_r_packages.R
# =============================================================================

BIOC_VERSION <- "3.23"   # pin: pairs with R 4.6.x

# ---- user-writable library (prepend, don't replace) -----------------------
user_lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(user_lib))
  stop("R_LIBS_USER is not set; cannot choose a user-writable library.")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
cat("Install target (.libPaths()[1]): ", .libPaths()[1], "\n", sep = "")

if (getRversion() < "4.6" || getRversion() >= "4.7")
  warning("This project pins Bioconductor ", BIOC_VERSION, ", which expects R 4.6.x; ",
          "you are on R ", as.character(getRversion()),
          ". Package versions may differ from the recorded environment.")

# ---- BiocManager, pinned --------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = user_lib, repos = "https://cloud.r-project.org")
if (as.character(BiocManager::version()) != BIOC_VERSION)
  BiocManager::install(version = BIOC_VERSION, ask = FALSE, update = FALSE)

# ---- helper: robust, non-interactive install ------------------------------
bioc_install <- function(pkgs)
  BiocManager::install(pkgs, lib = user_lib, update = FALSE, ask = FALSE,
                       INSTALL_opts = "--no-lock")

bioc_pkgs <- c(
  "minfi",
  "ChAMP",
  "limma",
  "DMRcate",                              # differentially methylated region calling
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "FlowSorted.Blood.EPIC",                # Houseman cell-type deconvolution reference
  "FlowSorted.CordBloodCombined.450k",    # dependency of FlowSorted.Blood.EPIC utils
  "sva",
  "clusterProfiler",
  "org.Hs.eg.db"
)

# Install per-package so one failure doesn't abort the rest; FlowSorted.Blood.EPIC
# (the historically fragile one) is retried and treated as fatal since the pipeline
# cannot do cell-type correction without it.
for (p in bioc_pkgs) {
  cat("\n---- installing:", p, "----\n")
  ok <- tryCatch({ if (!requireNamespace(p, quietly = TRUE)) bioc_install(p); TRUE },
                 error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); FALSE })
  if (!ok && p == "FlowSorted.Blood.EPIC") {
    cat("  retrying FlowSorted.Blood.EPIC...\n"); bioc_install(p)
  }
}

cran_pkgs <- c("tidyverse", "data.table", "matrixStats")
cran_need  <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(cran_need))
  install.packages(cran_need, lib = user_lib, repos = "https://cloud.r-project.org")

# ---- verification (fail loudly on a partial install) ----------------------
cat("\n=================== VERIFICATION ===================\n")
crit <- c("minfi", "limma", "DMRcate", "sva", "clusterProfiler",
          "IlluminaHumanMethylationEPICmanifest", "FlowSorted.Blood.EPIC")
missing <- crit[!vapply(crit, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing))
  stop("Critical packages NOT installed: ", paste(missing, collapse = ", "))

# The final run's Houseman step (01a_cell_counts.R) needs these exact objects:
suppressPackageStartupMessages(library(FlowSorted.Blood.EPIC))
needed <- c("projectCellType_CP", "IDOLOptimizedCpGs.compTable")
absent <- needed[!vapply(needed, exists, logical(1),
                         where = asNamespace("FlowSorted.Blood.EPIC"))]
if (length(absent))
  stop("FlowSorted.Blood.EPIC is installed but missing objects the pipeline uses: ",
       paste(absent, collapse = ", "),
       " -- 01a_cell_counts.R would fail. Check the package version.")
cat("FlowSorted.Blood.EPIC objects for 01a present:", paste(needed, collapse = ", "), "\n")

# ---- record the environment -----------------------------------------------
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)
manifest <- file.path("results/qc", "r_package_versions.txt")
vers <- c(paste("Recorded by install_r_packages.R"),
          paste("R:", as.character(getRversion())),
          paste("Bioconductor:", as.character(BiocManager::version())),
          "", "Critical package versions:")
for (p in unique(c(crit, "FlowSorted.CordBloodCombined.450k",
                   "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "matrixStats")))
  vers <- c(vers, sprintf("  %-46s %s", p,
                          tryCatch(as.character(packageVersion(p)),
                                   error = function(e) "NOT INSTALLED")))
writeLines(vers, manifest)
cat("\nWrote version manifest -> ", manifest, "\n", sep = "")
cat("Setup complete and verified.\n")
