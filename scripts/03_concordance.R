# =============================================================================
# 03_concordance.R
#
# Diet-quality methylation<->expression concordance analysis (LIFE-Adult, n=48).
#
# RESEARCH QUESTION
#   Which CpG methylation changes associated with diet quality (continuous
#   DietScore) show CONCORDANT changes in expression of the corresponding gene,
#   independent of BMI, age, and sex?
#
# THREE-SIGNAL CONCORDANCE FRAMEWORK
#   (1) Diet -> methylation   EWAS:  Mval ~ DietScore + covars + cell-type
#   (2) Diet -> expression    DE:    expr ~ DietScore + covars
#   (3) Methylation <-> expr  eQTM:  expr ~ Mval + covars   (local link)
#   A CpG-gene pair is CONCORDANT when the diet->expression effect matches what
#   the diet->methylation effect predicts THROUGH the local eQTM slope:
#       sign(b_diet_expr)  ==  sign(b_diet_meth * b_eqtm)
#   This is data-driven (does not assume promoter=inverse / body=positive) but
#   the CpG's gene region is retained for interpretation.
#
# POWER NOTE (n=48)
#   Genome-wide FDR across ~750k probes is underpowered -> EWAS is treated as
#   DISCOVERY (nominal-p threshold + top-N), and multiple testing is controlled
#   WITHIN the resulting CpG-gene pair set. A candidate-gene mode (cfg$candidate_*)
#   supports a properly-powered, a-priori-restricted test in parallel.
#
# DESIGN NOTE -- eQTM p-value is a GATE, not an FDR-tested hypothesis (deliberate).
#   The BH-FDR (cfg$concordance_fdr) is applied to the two DIET signals per pair
#   (pmax of diet->meth and diet->expr p). The local methylation<->expression link
#   (eqtm_p) is used only as a hard filter on which pairs qualify as concordant
#   (cfg$eqtm_p_max). It is intentionally kept OUT of the FDR set: it is a
#   supporting/mechanistic condition, not one of the diet-association hypotheses we
#   are counting. Folding it in would conflate "diet moves this" with "methylation
#   locally tracks expression" and distort the correction. This is a choice, not an
#   oversight -- revisit if a reviewer wants a single combined error rate.
#
# Inputs  (from scripts 01 & 02):
#   data/processed/m_values.rds, pheno_methylation.rds   (has 6 cell-type props)
#   data/processed/grSet_normalized_filtered.rds         (for CpG annotation)
#   data/processed/expr_gene_paired.rds, pheno_expression.rds
# Outputs : results/tables/*.csv, results/figures/*.pdf
#
# STATUS: SKETCH. The methylation side is fully specified; the expression side
#         depends on script-02 outputs that require the not-yet-downloaded
#         values file. Run once 01 & 02 have produced their processed objects.
#
# Run non-interactively:  Rscript scripts/03_concordance.R
# =============================================================================

suppressPackageStartupMessages({
  library(minfi)                                       # getAnnotation on grSet
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(limma)
  library(data.table)
})

# ----------------------------------------------------------------------------
# 0. Configuration
# ----------------------------------------------------------------------------
cfg <- list(
  proc_dir   = "data/processed",
  tbl_dir    = "results/tables",
  fig_dir    = "results/figures",

  # analysis variable + adjustment set
  exposure   = "DietScore",
  covariates = c("Age", "Sex", "BMI"),
  # 6 Houseman proportions from script 01; ONE dropped (they sum to ~1 -> collinear
  # with the intercept). Neu (neutrophils, most abundant) is the dropped reference.
  cell_types      = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
  cell_ref_drop   = "Neu",
  adjust_expr_for_cells = FALSE,   # blood expr: usually adjust clinical covars only

  # discovery thresholds
  ewas_p_discovery = 1e-3,   # nominal p to carry a CpG into the concordance stage
  ewas_top_n       = 5000,   # ...or top-N by p, whichever is larger (union)
  concordance_fdr  = 0.05,   # BH-FDR applied WITHIN the CpG-gene pair set
  eqtm_p_max       = 0.05,   # local methylation-expression link must be nominal

  # candidate-gene mode (optional a-priori list; one gene symbol per line)
  candidate_genes_file = "data/raw/annotation/candidate_genes.txt",

  # CpG annotation columns (EPIC ilm10b4.hg19)
  gene_col   = "UCSC_RefGene_Name",   # semicolon-separated gene symbols
  group_col  = "UCSC_RefGene_Group",  # TSS1500;Body;5'UTR;...

  seed = 1234
)

set.seed(cfg$seed)
dir.create(cfg$tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$fig_dir, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))
rds <- function(f) { p <- file.path(cfg$proc_dir, f)
  if (!file.exists(p)) stop("Missing input: ", p, " -- run scripts 01/02 first.")
  readRDS(p) }

# ----------------------------------------------------------------------------
# 1. Load processed data & harmonize to the shared sample set
# ----------------------------------------------------------------------------
mval       <- rds("m_values.rds")                 # CpG x sample (M-values)
pheno_meth <- rds("pheno_methylation.rds")         # has DietScore, Age, Sex, BMI, cells
grSet      <- rds("grSet_normalized_filtered.rds") # for CpG->gene annotation
expr_gene  <- rds("expr_gene_paired.rds")          # gene x sample (log2), cols = Sample_Name

# Common samples (both assays already key on methylation Sample_Name).
common <- Reduce(intersect, list(colnames(mval), colnames(expr_gene),
                                 rownames(pheno_meth)))
log_msg("Samples: meth=%d, expr=%d -> analysing %d paired.",
        ncol(mval), ncol(expr_gene), length(common))
if (length(common) < 10) stop("Too few paired samples (", length(common), ").")

mval      <- mval[, common, drop = FALSE]
expr_gene <- expr_gene[, common, drop = FALSE]
ph        <- pheno_meth[common, , drop = FALSE]
ph$Sex    <- droplevels(factor(ph$Sex))

# ----------------------------------------------------------------------------
# 2. Design matrices
# ----------------------------------------------------------------------------
cells_used <- setdiff(cfg$cell_types, cfg$cell_ref_drop)
build_design <- function(ph, with_cells) {
  terms <- c(cfg$exposure, cfg$covariates, if (with_cells) cells_used)
  terms <- terms[terms %in% names(ph)]
  f <- as.formula(paste("~", paste(terms, collapse = " + ")))
  mm <- model.matrix(f, data = ph)
  list(mm = mm, coef = make.names(cfg$exposure))
}
d_meth <- build_design(ph, with_cells = TRUE)
d_expr <- build_design(ph, with_cells = cfg$adjust_expr_for_cells)
log_msg("Methylation model:  ~ %s", paste(colnames(d_meth$mm)[-1], collapse = " + "))
log_msg("Expression  model:  ~ %s", paste(colnames(d_expr$mm)[-1], collapse = " + "))
if (nrow(d_meth$mm) < length(common))
  log_msg("WARNING: %d samples dropped by model.matrix (missing covariates).",
          length(common) - nrow(d_meth$mm))

# ----------------------------------------------------------------------------
# 3. (1) Diet -> methylation  EWAS  (limma, moderated t)
# ----------------------------------------------------------------------------
log_msg("EWAS: %s ~ diet + covars + cell-type over %d CpGs...",
        "M-value", nrow(mval))
fit_m  <- eBayes(lmFit(mval[, rownames(d_meth$mm)], d_meth$mm))
ewas   <- topTable(fit_m, coef = d_meth$coef, number = Inf, sort.by = "p")
ewas$CpG <- rownames(ewas)
# Full table: gzipped (~35 MB), GITIGNORED + regenerable. fread() auto-decompresses.
# (03a_ewas.R writes this same pair of files for EWAS-only runs.)
fwrite(ewas, file.path(cfg$tbl_dir, "ewas_diet_methylation.csv.gz"))
# Tracked artifact: the nominal-significant discovery pool only (small).
sig <- ewas[ewas$P.Value < cfg$ewas_p_discovery, ]
fwrite(sig[order(sig$P.Value), ], file.path(cfg$tbl_dir, "ewas_diet_significant.csv"))
log_msg("EWAS done. Genome-wide FDR<0.05: %d CpGs (expected ~0 at n=48).",
        sum(ewas$adj.P.Val < 0.05))

# ----------------------------------------------------------------------------
# 4. CpG -> gene mapping (EPIC annotation), restricted to measured genes
# ----------------------------------------------------------------------------
ann <- as.data.frame(getAnnotation(grSet))[, c(cfg$gene_col, cfg$group_col)]
ann$CpG <- rownames(ann)
# explode semicolon-separated gene / region into one row per CpG-gene pair
explode <- function(df) {
  g  <- strsplit(df[[cfg$gene_col]],  ";", fixed = TRUE)
  gr <- strsplit(df[[cfg$group_col]], ";", fixed = TRUE)
  n  <- lengths(g)
  data.frame(CpG = rep(df$CpG, n),
             gene   = unlist(g),
             region = unlist(Map(function(a, k) if (length(a)) a else rep(NA, k), gr, n)),
             stringsAsFactors = FALSE)
}
cpg2gene <- unique(explode(ann[ann[[cfg$gene_col]] != "", ]))
cpg2gene <- cpg2gene[cpg2gene$gene %in% rownames(expr_gene), ]
log_msg("CpG-gene pairs with a measured gene: %d (%d unique genes).",
        nrow(cpg2gene), length(unique(cpg2gene$gene)))

# ----------------------------------------------------------------------------
# 5. Select the discovery / candidate pair set
# ----------------------------------------------------------------------------
disc_cpgs <- union(
  ewas$CpG[ewas$P.Value < cfg$ewas_p_discovery],
  head(ewas$CpG[order(ewas$P.Value)], cfg$ewas_top_n)
)
candidate_genes <- character(0)
if (file.exists(cfg$candidate_genes_file)) {
  raw <- readLines(cfg$candidate_genes_file)
  raw <- sub("#.*$", "", raw)                    # strip full-line and inline comments
  candidate_genes <- unique(trimws(raw))
  candidate_genes <- candidate_genes[candidate_genes != ""]
  log_msg("Candidate-gene list: %d genes.", length(candidate_genes))
}
pairs <- cpg2gene[cpg2gene$CpG %in% disc_cpgs |
                  cpg2gene$gene %in% candidate_genes, ]
pairs$candidate <- pairs$gene %in% candidate_genes
log_msg("Concordance pair set: %d pairs (%d from candidate genes).",
        nrow(pairs), sum(pairs$candidate))
if (!nrow(pairs)) stop("No CpG-gene pairs selected -- loosen ewas_p_discovery.")

# ----------------------------------------------------------------------------
# 6. (2) Diet -> expression  DE   for the genes in the pair set
# ----------------------------------------------------------------------------
genes_test <- unique(pairs$gene)
fit_e <- eBayes(lmFit(expr_gene[genes_test, rownames(d_expr$mm), drop = FALSE], d_expr$mm))
de <- topTable(fit_e, coef = d_expr$coef, number = Inf, sort.by = "none")
de$gene <- rownames(de)
fwrite(de, file.path(cfg$tbl_dir, "de_diet_expression.csv"))
log_msg("DE done for %d genes.", nrow(de))

# ----------------------------------------------------------------------------
# 7. (3) eQTM  +  concordance integration, per CpG-gene pair
# ----------------------------------------------------------------------------
# eQTM slope: expr_gene[gene] ~ Mval[CpG] + covariates, over the paired samples.
cov_mm <- d_expr$mm[, setdiff(colnames(d_expr$mm), c("(Intercept)", d_expr$coef)),
                    drop = FALSE]
samp   <- rownames(d_expr$mm)
na_eqtm <- c(slope = NA_real_, p = NA_real_)
eqtm_slope <- function(cpg, gene) {
  tryCatch({
    y <- as.numeric(expr_gene[gene, samp])
    x <- as.numeric(mval[cpg, samp])
    X <- cbind(1, x, cov_mm)
    # not estimable: missing values, no methylation variance, or rank-deficient
    if (anyNA(y) || anyNA(x) || stats::sd(x) == 0) return(na_eqtm)
    f <- lm.fit(X, y)                            # coef 2 = methylation slope
    if (f$rank < ncol(X)) return(na_eqtm)
    b   <- f$coefficients[2]
    rss <- sum(f$residuals^2); df <- length(y) - f$rank
    se  <- sqrt(diag(chol2inv(f$qr$qr[1:f$rank, 1:f$rank, drop = FALSE])) * rss / df)[2]
    c(slope = b, p = 2 * pt(-abs(b / se), df))
  }, error = function(e) na_eqtm)
}
log_msg("Computing eQTM for %d pairs...", nrow(pairs))
eq <- t(mapply(eqtm_slope, pairs$CpG, pairs$gene))
pairs$eqtm_slope <- eq[, "slope"]
pairs$eqtm_p     <- eq[, "p"]
n_bad <- sum(is.na(pairs$eqtm_p))
if (n_bad) log_msg("  %d/%d pairs non-estimable (eQTM = NA), treated as non-concordant.",
                   n_bad, nrow(pairs))

# attach diet effects
mi <- match(pairs$CpG, ewas$CpG)
pairs$diet_meth_b <- ewas$logFC[mi];  pairs$diet_meth_p <- ewas$P.Value[mi]
gi <- match(pairs$gene, de$gene)
pairs$diet_expr_b <- de$logFC[gi];    pairs$diet_expr_p <- de$P.Value[gi]

# concordance: does the diet->expr effect match diet->meth propagated via eQTM?
pairs$predicted_expr_sign <- sign(pairs$diet_meth_b * pairs$eqtm_slope)
pairs$concordant <- with(pairs,
  sign(diet_expr_b) == predicted_expr_sign &
  eqtm_p     <= cfg$eqtm_p_max &
  diet_meth_p <= cfg$ewas_p_discovery)
pairs$concordant[is.na(pairs$concordant)] <- FALSE   # non-estimable eQTM -> not concordant

# BH-FDR within the pair set, on the weakest of the two diet signals per pair
pairs$pair_p   <- pmax(pairs$diet_meth_p, pairs$diet_expr_p)
pairs$pair_fdr <- p.adjust(pairs$pair_p, method = "BH")
pairs$hit      <- pairs$concordant & pairs$pair_fdr <= cfg$concordance_fdr

pairs <- pairs[order(!pairs$hit, pairs$pair_p), ]
fwrite(pairs, file.path(cfg$tbl_dir, "concordance_pairs.csv"))
log_msg("Concordant pairs: %d nominal, %d at FDR<%.2f (%d in candidate genes).",
        sum(pairs$concordant), sum(pairs$hit), cfg$concordance_fdr,
        sum(pairs$hit & pairs$candidate))

# ----------------------------------------------------------------------------
# 8. Figures
# ----------------------------------------------------------------------------
# EWAS volcano (rasterized PNG: 806k points -> KB not a ~31 MB vector PDF)
png(file.path(cfg$fig_dir, "ewas_volcano.png"), width = 1200, height = 1200, res = 200)
plot(ewas$logFC, -log10(ewas$P.Value), pch = 16, cex = 0.3, col = "grey60",
     xlab = "Diet effect on M-value", ylab = "-log10 p", main = "EWAS: diet -> methylation")
abline(h = -log10(cfg$ewas_p_discovery), col = "red", lty = 2)
dev.off()

# Concordance scatter: diet->meth effect vs diet->expr effect for pair set
pdf(file.path(cfg$fig_dir, "concordance_scatter.pdf"), width = 6, height = 6)
col <- ifelse(pairs$hit, "firebrick", ifelse(pairs$concordant, "orange", "grey70"))
plot(pairs$diet_meth_b, pairs$diet_expr_b, pch = 16, cex = 0.5, col = col,
     xlab = "Diet -> methylation effect", ylab = "Diet -> expression effect",
     main = "Concordance (red = FDR hit)")
abline(h = 0, v = 0, col = "grey85"); dev.off()

writeLines(capture.output(sessionInfo()),
           file.path(cfg$tbl_dir, "sessionInfo_03_concordance.txt"))
log_msg("Done. Tables in %s/, figures in %s/.", cfg$tbl_dir, cfg$fig_dir)
