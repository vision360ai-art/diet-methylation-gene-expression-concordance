#!/usr/bin/env Rscript
# One-off: confirm data/processed/m_values.rds is a valid, complete M-value matrix.
p <- "data/processed/m_values.rds"
cat("file exists:", file.exists(p), "\n")
cat("size_MB:", round(file.info(p)$size / 1024^2, 1), "\n")
m <- readRDS(p)
cat("class:", class(m)[1], "| storage:", typeof(m), "\n")
cat("dim:", paste(dim(m), collapse = " x "), "\n")
cat("anyNA:", anyNA(m), "| anyInf:", any(is.infinite(m)), "\n")
cat("range:", paste(round(range(m), 3), collapse = " .. "), "\n")
cat("rownames[1:2]:", paste(head(rownames(m), 2), collapse = ", "), "\n")
cat("colnames[1:3]:", paste(head(colnames(m), 3), collapse = ", "), "\n")
cat("OK\n")
