#!/usr/bin/env Rscript
# Install ONLY FlowSorted.Blood.EPIC into the USER-WRITABLE R library.
#
# This is a targeted helper -- the canonical, fully-verified setup is
# scripts/install_r_packages.R (it installs everything, pins Bioc, and checks the
# exact objects the pipeline uses). Use this when only the cell-type reference needs
# (re)installing. The two MUST resolve the same version, so the Bioc pin is shared.
#
# The default .libPaths() first entry is C:/Program Files/R/.../library, which
# is not writable without admin rights. We target R_LIBS_USER instead — same fix
# used for BiocManager earlier. Run: Rscript scripts/install_flowsorted.R

BIOC_VERSION <- "3.23"   # keep in sync with install_r_packages.R

user_lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(user_lib)) stop("R_LIBS_USER is not set; cannot locate user library.")

dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))          # user lib FIRST so installs land there

cat("Install target (.libPaths()[1]):", .libPaths()[1], "\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib, repos = "https://cloud.r-project.org")
}
if (as.character(BiocManager::version()) != BIOC_VERSION)
  BiocManager::install(version = BIOC_VERSION, ask = FALSE, update = FALSE)

# FlowSorted.Blood.EPIC provides the IDOL reference. NOTE: the pipeline does NOT
# call estimateCellCounts2 (it co-normalizes with the full reference and OOMs on
# this host); 01a_cell_counts.R uses the lighter projectCellType_CP +
# IDOLOptimizedCpGs.compTable path instead. So we verify THOSE objects below.
BiocManager::install(
  "FlowSorted.Blood.EPIC",
  lib          = user_lib,
  update       = FALSE,
  ask          = FALSE,
  INSTALL_opts = "--no-lock"
)

# ---- confirm it loads from the user library ------------------------------
ok <- requireNamespace("FlowSorted.Blood.EPIC", quietly = TRUE)
cat("\nFlowSorted.Blood.EPIC installed & loadable:", ok, "\n")
if (!ok) stop("Install did not produce a loadable package.")
cat("Loaded from:", dirname(dirname(getNamespaceInfo("FlowSorted.Blood.EPIC", "path"))), "\n")
cat("Version:", as.character(packageVersion("FlowSorted.Blood.EPIC")), "\n")

needed <- c("projectCellType_CP", "IDOLOptimizedCpGs.compTable")
absent <- needed[!vapply(needed, exists, logical(1),
                         where = asNamespace("FlowSorted.Blood.EPIC"))]
if (length(absent))
  stop("Installed, but 01a's required objects are missing: ",
       paste(absent, collapse = ", "))
cat("Objects required by 01a_cell_counts.R present:", paste(needed, collapse = ", "), "\n")
