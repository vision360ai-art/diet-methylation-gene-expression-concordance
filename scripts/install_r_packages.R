# R environment setup for diet-methylation-expression project
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
.libPaths(Sys.getenv("R_LIBS_USER"))

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "minfi",
  "ChAMP",
  "limma",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "sva",
  "clusterProfiler",
  "org.Hs.eg.db"
))

install.packages(c("tidyverse", "data.table"))

cat("Setup complete. Run sessionInfo() to confirm package versions.\n")
