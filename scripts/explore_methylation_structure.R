#!/usr/bin/env Rscript
# Exploratory visualization of the 36 preprocessed methylation samples, BEFORE and
# AFTER adjusting for blood cell composition (the dominant whole-blood confounder).
#
#   RAW      : PCA / clustering on the most-variable CpGs of the normalized beta.
#   CORRECTED: same, after residualizing every CpG on the 6 Houseman cell-type
#              proportions (5 used as covariates; one dropped as reference because
#              proportions sum to ~1). This is the "cell-corrected" view -- it strips
#              the leukocyte-composition signal so any remaining structure is not
#              just cell-mix.
#
# Note: the beta matrix is identical to the pre-Houseman run (normalization does not
# use cell counts) -- the correction is applied HERE by regressing the proportions
# out, using the cell-type columns now present in pheno.
# Exploratory only; not part of the numbered pipeline.

suppressMessages({
  library(ggplot2)
  library(matrixStats)
})

cfg <- list(
  beta      = "data/processed/beta_values.rds",
  pheno     = "data/processed/pheno_methylation.csv",
  outdir    = "results/figures",
  n_top     = 20000L,   # most-variable CpGs used for PCA/clustering
  cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
)
dir.create(cfg$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- load ----------------------------------------------------------------
beta  <- readRDS(cfg$beta)                       # CpG x sample
pheno <- read.csv(cfg$pheno, check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(all(colnames(beta) %in% pheno$Sample_Name))
pheno <- pheno[match(colnames(beta), pheno$Sample_Name), ]
cat(sprintf("Loaded beta: %d CpGs x %d samples\n", nrow(beta), ncol(beta)))

have_ct <- all(cfg$cellTypes %in% names(pheno))
if (!have_ct)
  stop("pheno lacks cell-type columns (", paste(cfg$cellTypes, collapse=","),
       ") -- run 01a_cell_counts.R + 01 so the correction is possible.")
cat("Cell-type proportions present. Means: ",
    paste(sprintf("%s=%.3f", cfg$cellTypes, colMeans(pheno[, cfg$cellTypes])),
          collapse = ", "), "\n", sep = "")

# ---- cell-composition residualization ------------------------------------
# Residualize each CpG on 5 cell-type proportions (drop 1 as reference to avoid the
# sum-to-one collinearity). Projection is a cheap 36x36 hat-matrix applied to the
# beta matrix: R = Y %*% (I - H), H = X (X'X)^-1 X'. Subtracting a per-CpG constant
# across samples leaves PCA (which centers) and sample-distance (shift-invariant)
# unchanged, so removing the intercept term with the fit is harmless here.
X   <- cbind(Intercept = 1, as.matrix(pheno[, cfg$cellTypes[-length(cfg$cellTypes)]]))
H   <- X %*% solve(crossprod(X)) %*% t(X)          # 36 x 36
I_H <- diag(ncol(beta)) - H
beta_corr <- beta %*% I_H                          # CpG x sample, cell-adjusted
colnames(beta_corr) <- colnames(beta)
cat("Residualized beta on 5 cell-type proportions.\n")

# ---- helper: top-variable subset + PCA -----------------------------------
top_pca <- function(mat, n_top) {
  v   <- matrixStats::rowVars(mat)
  top <- order(v, decreasing = TRUE)[seq_len(min(n_top, nrow(mat)))]
  matT <- mat[top, ]
  pca  <- prcomp(t(matT), center = TRUE, scale. = FALSE)
  ve   <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  list(matT = matT, pca = pca, ve = ve, ntop = length(top))
}

make_pcdf <- function(res) data.frame(
  Sample_Name  = rownames(res$pca$x),
  PC1 = res$pca$x[, 1], PC2 = res$pca$x[, 2], PC3 = res$pca$x[, 3],
  DietScore    = pheno$DietScore,
  HLS_group    = pheno$HLS_group,
  Sample_Group = pheno$Sample_Group,
  Sex          = pheno$Sex,
  Age          = pheno$Age,
  BMI          = pheno$BMI,
  Neu          = pheno$Neu,
  CD4T         = pheno$CD4T
)

raw  <- top_pca(beta,      cfg$n_top)
corr <- top_pca(beta_corr, cfg$n_top)
pdf_raw  <- make_pcdf(raw)
pdf_corr <- make_pcdf(corr)

# ---- plotting helpers -----------------------------------------------------
lab <- function(ve, i) sprintf("PC%d (%.1f%%)", i, ve[i])

pca_cont <- function(d, ve, colvar, title, sub, opt = "C") {
  ggplot(d, aes(PC1, PC2, color = .data[[colvar]])) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_viridis_c(option = opt) +
    labs(title = title, subtitle = sub, x = lab(ve, 1), y = lab(ve, 2),
         color = colvar) +
    theme_bw(base_size = 12)
}
pca_cat <- function(d, ve, colvar, title, sub, cols = NULL) {
  g <- ggplot(d, aes(PC1, PC2, color = .data[[colvar]])) +
    geom_point(size = 3, alpha = 0.9) +
    labs(title = title, subtitle = sub, x = lab(ve, 1), y = lab(ve, 2)) +
    theme_bw(base_size = 12)
  if (!is.null(cols)) g <- g + scale_color_manual(values = cols)
  g
}

hls_cols <- c(low = "#D55E00", high = "#0072B2")
sub_raw  <- sprintf("RAW beta | top %d variable CpGs, n=%d", raw$ntop,  ncol(beta))
sub_corr <- sprintf("CELL-CORRECTED beta | top %d variable CpGs, n=%d", corr$ntop, ncol(beta))

# RAW: colored by neutrophil proportion -> shows cell composition IS a top driver
p1 <- pca_cont(pdf_raw, raw$ve, "Neu",
               "Methylation PCA (RAW) — colored by neutrophil proportion",
               sub_raw, opt = "D")
# RAW: by DietScore / HLS (for before/after comparison)
p2 <- pca_cont(pdf_raw, raw$ve, "DietScore",
               "Methylation PCA (RAW) — colored by DietScore", sub_raw)
p3 <- pca_cat(pdf_raw, raw$ve, "HLS_group",
              "Methylation PCA (RAW) — colored by lifestyle (HLS) group",
              sub_raw, hls_cols)

# CORRECTED: by neutrophil (should now be washed out), DietScore, HLS
p4 <- pca_cont(pdf_corr, corr$ve, "Neu",
               "Methylation PCA (CELL-CORRECTED) — colored by neutrophil proportion",
               sub_corr, opt = "D")
p5 <- pca_cont(pdf_corr, corr$ve, "DietScore",
               "Methylation PCA (CELL-CORRECTED) — colored by DietScore", sub_corr)
p6 <- pca_cat(pdf_corr, corr$ve, "HLS_group",
              "Methylation PCA (CELL-CORRECTED) — colored by lifestyle (HLS) group",
              sub_corr, hls_cols)

ggsave(file.path(cfg$outdir, "pca_raw_by_neutrophil.png"),        p1, width = 7, height = 5.5, dpi = 150)
ggsave(file.path(cfg$outdir, "pca_raw_by_dietscore.png"),         p2, width = 7, height = 5.5, dpi = 150)
ggsave(file.path(cfg$outdir, "pca_raw_by_hlsgroup.png"),          p3, width = 7, height = 5.5, dpi = 150)
ggsave(file.path(cfg$outdir, "pca_cellcorr_by_neutrophil.png"),   p4, width = 7, height = 5.5, dpi = 150)
ggsave(file.path(cfg$outdir, "pca_cellcorr_by_dietscore.png"),    p5, width = 7, height = 5.5, dpi = 150)
ggsave(file.path(cfg$outdir, "pca_cellcorr_by_hlsgroup.png"),     p6, width = 7, height = 5.5, dpi = 150)

# ---- hierarchical clustering (cell-corrected) ----------------------------
dend <- function(res, fname, title) {
  d  <- dist(t(res$matT)); hc <- hclust(d, method = "ward.D2")
  png(file.path(cfg$outdir, fname), width = 1400, height = 800, res = 150)
  op <- par(mar = c(7, 4, 3, 1))
  lbl <- sprintf("%s [%s/D%s]", pheno$Sample_Name, pheno$HLS_group, pheno$DietScore)
  plot(hc, labels = lbl, hang = -1, cex = 0.7, main = title, xlab = "", sub = "")
  par(op); dev.off()
}
dend(raw,  "hclust_raw_dendrogram.png",
     "Hierarchical clustering — RAW beta (Ward.D2, top variable CpGs)")
dend(corr, "hclust_cellcorr_dendrogram.png",
     "Hierarchical clustering — CELL-CORRECTED beta (Ward.D2, top variable CpGs)")

# ---- before/after association tables -------------------------------------
assoc <- function(d, tag) {
  cat(sprintf("\n=== %s: PC correlation with continuous vars ===\n", tag))
  num <- d[, c("PC1","PC2","PC3","DietScore","Age","BMI","Neu","CD4T")]
  print(round(cor(num, use = "pairwise.complete.obs")[
      c("DietScore","Age","BMI","Neu","CD4T"), c("PC1","PC2","PC3")], 3))
  r2 <- function(pc, grp) summary(lm(pc ~ grp))$r.squared
  cat(sprintf("%s: categorical variance explained (R^2 of PC ~ group):\n", tag))
  for (g in c("HLS_group","Sample_Group","Sex")) {
    cat(sprintf("  %-13s  PC1=%.3f  PC2=%.3f  PC3=%.3f\n", g,
                r2(d$PC1, d[[g]]), r2(d$PC2, d[[g]]), r2(d$PC3, d[[g]])))
  }
}
assoc(pdf_raw,  "RAW")
assoc(pdf_corr, "CELL-CORRECTED")

cat(sprintf("\nFigures written to %s/\n", cfg$outdir))
