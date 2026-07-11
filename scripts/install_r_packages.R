# R environment setup for diet-methylation-expression project
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
.libPaths(Sys.getenv("R_LIBS_USER"))

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "minfi",
  "ChAMP",
  "limma",
  "DMRcate",                       # differentially methylated region calling
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "FlowSorted.Blood.EPIC",         # Houseman cell-type deconvolution reference
  "FlowSorted.CordBloodCombined.450k",  # dependency of FlowSorted.Blood.EPIC utils
  "sva",
  "clusterProfiler",
  "org.Hs.eg.db"
))

install.packages(c("tidyverse", "data.table"))

cat("Setup complete. Run sessionInfo() to confirm package versions.\n")
