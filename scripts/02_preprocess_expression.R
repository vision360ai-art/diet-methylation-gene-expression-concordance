# =============================================================================
# 02_preprocess_expression.R
#
# Preprocessing of Illumina HumanHT-12 v4 expression data for the LIFE-Adult
# diet / methylation / expression concordance project.
#
# CONTEXT
#   The distributed LIFE expression dataset is delivered ALREADY normalized,
#   batch-corrected and imputed (see Gx_DataDictionary: "imputation of this
#   probe in all samples", "even after batch correction"). This script therefore
#   LOADS a processed probe x sample matrix and applies the provider's per-probe
#   QC flags -- it does NOT re-run bead-level background correction / quantile
#   normalization. (If you instead receive a raw GenomeStudio "Sample Probe
#   Profile", that needs a different pipeline: limma::neqc -> filter -> ComBat.)
#
#   Sample linkage across assays:
#       expression `sampleID`  ==  methylation `Synonym`  (SampleSheet_V2.csv)
#   Confirmed overlap = 48 subjects. Output is scoped to those 48 paired samples,
#   relabelled to the methylation `Sample_Name` so both assays share one key.
#
# Pipeline:
#   1. Build the cross-assay sample linkage table (n=48 paired)
#   2. Load the processed expression matrix
#   3. Load the probe annotation table (gene mapping + QC flags)
#   4. log2 sanity check
#   5. Probe-level QC filtering (provider flags; each self-skips if absent)
#   6. Subset to n=48 paired samples, relabel columns -> Sample_Name
#   7. (optional) ComBat batch correction  -- OFF by default (already corrected)
#   8. (optional) gene-level collapse
#   9. QC figures (density, PCA by batch) + save outputs
#
# Outputs : data/processed/expr_*.rds, data/processed/pheno_expression.*,
#           results/qc/qc_expr_*.pdf
#
# Run non-interactively:  Rscript scripts/02_preprocess_expression.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)   # ComBat lives in sva; limma used for plotMDS / general utils
})

# ----------------------------------------------------------------------------
# 0. Configuration
#
#  >>> CONFIRM the marked (# CONFIRM) values once the expression VALUES file is
#      downloaded -- filename and column names are not yet known. Everything
#      else follows the Gx_DataDictionary spec. <<<
# ----------------------------------------------------------------------------
cfg <- list(
  # --- inputs ---
  sample_sheet   = "data/raw/SampleSheet_V2.csv",              # methylation sheet (linkage)
  sample_info    = "data/raw/Gx_sample_info-21-06-15_v3.txt",  # expression sample metadata
  expr_matrix    = "data/raw/Gx_expression_matrix.txt",        # CONFIRM: probe x sample values (STILL MISSING)
  probe_annot    = "data/raw/Gx_probe_info-21-06-15_v3.txt",   # confirmed: 48106 probes x 42 cols

  out_dir        = "data/processed",
  qc_dir         = "results/qc",

  # --- column names ---
  link_meth_key  = "Synonym",       # column in SampleSheet_V2.csv that == expression sampleID
  link_meth_name = "Sample_Name",   # methylation primary sample key (target label)
  expr_id_col    = "sampleID",      # sample-id column in sample_info
  probe_id_col   = "PROBE_ID",      # CONFIRM: probe-id column in expr_matrix & probe_annot
  gene_symbol_col= "symbol_INGENUITY",
  gene_entrez_col= "ilmn_entrezID_INGENUITY",

  # --- probe QC flags (from probe_annot; each filter self-skips if col absent) ---
  flag_qcok      = "probe_QCok",              # keep TRUE
  flag_control   = "is_purecontrol",          # drop TRUE
  flag_expressed = "expressed_blood",         # keep TRUE
  flag_bestprobe = "bestexpressed_probe_in_gene",  # used only for gene collapse
  batch_assoc_cols = c("strong_processing_batch_association",
                       "strong_temporal_batch_association",
                       "strong_hybridisierungchip_batch_association"),

  # --- switches ---
  probes_in_rows = TRUE,     # TRUE: matrix is probes x samples; FALSE: transpose on load
  already_log2   = NA,       # NA = auto-detect from value range; TRUE/FALSE to force
  run_combat     = FALSE,    # data ships batch-corrected; leave OFF unless QC says otherwise
  combat_batch   = "processing_batch",   # batch var if run_combat = TRUE
  collapse_genes = TRUE,     # also emit a gene-level matrix
  best_sample_only = TRUE,   # keep only BEST_SAMPLE == TRUE rows of sample_info

  seed = 1234
)

set.seed(cfg$seed)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$qc_dir,  recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))

# Robust coercion of provider TRUE/FALSE flags (may be logical, "TRUE", "1", "yes").
as_flag <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(trimws(as.character(x)))
  out <- x %in% c("true", "t", "1", "yes", "y")
  out[is.na(x) | x == ""] <- NA
  out
}

need_file <- function(path, what) {
  if (!file.exists(path))
    stop(sprintf("%s not found: '%s'. Download it or fix the cfg path.", what, path))
}

# ----------------------------------------------------------------------------
# 1. Cross-assay sample linkage table (n=48 paired)
# ----------------------------------------------------------------------------
need_file(cfg$sample_sheet, "Methylation sample sheet")
need_file(cfg$sample_info,  "Expression sample_info")

meth <- read.csv(cfg$sample_sheet, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c(cfg$link_meth_key, cfg$link_meth_name) %in% names(meth)))
  stop("SampleSheet missing linkage columns: ",
       paste(setdiff(c(cfg$link_meth_key, cfg$link_meth_name), names(meth)), collapse = ", "))

sinfo <- fread(cfg$sample_info, data.table = FALSE)
if (cfg$best_sample_only && "BEST_SAMPLE" %in% names(sinfo)) {
  keep <- as_flag(sinfo$BEST_SAMPLE)
  log_msg("sample_info: keeping %d/%d BEST_SAMPLE==TRUE", sum(keep, na.rm = TRUE), nrow(sinfo))
  sinfo <- sinfo[which(keep), ]
}

# link: expression sampleID == methylation Synonym
link <- merge(
  meth[, c(cfg$link_meth_name, cfg$link_meth_key)],
  sinfo,
  by.x = cfg$link_meth_key, by.y = cfg$expr_id_col,
  all = FALSE
)
# recover the expression id under a stable name
link[[cfg$expr_id_col]] <- link[[cfg$link_meth_key]]
link <- link[!duplicated(link[[cfg$expr_id_col]]), ]
log_msg("Paired samples (methylation Synonym == expression sampleID): n = %d",
        nrow(link))
if (nrow(link) == 0) stop("No paired samples found -- check linkage columns.")

# ----------------------------------------------------------------------------
# 2. Load processed expression matrix
# ----------------------------------------------------------------------------
need_file(cfg$expr_matrix, "Expression values matrix")
log_msg("Loading expression matrix: %s", cfg$expr_matrix)
em <- fread(cfg$expr_matrix, data.table = FALSE)

if (!cfg$probes_in_rows) {           # samples x probes -> transpose to probes x samples
  rn <- em[[1]]; em <- as.data.frame(t(em[, -1])); colnames(em) <- rn
  em <- cbind(data.frame(probe = rownames(em)), em); names(em)[1] <- cfg$probe_id_col
}
if (!cfg$probe_id_col %in% names(em))
  stop(sprintf("Probe-id column '%s' not in expression matrix. Got: %s",
               cfg$probe_id_col, paste(head(names(em)), collapse = ", ")))

probe_ids <- as.character(em[[cfg$probe_id_col]])
expr <- as.matrix(em[, setdiff(names(em), cfg$probe_id_col), drop = FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- probe_ids
log_msg("Expression matrix: %d probes x %d samples", nrow(expr), ncol(expr))

# ----------------------------------------------------------------------------
# 3. Load probe annotation (gene mapping + QC flags)
# ----------------------------------------------------------------------------
need_file(cfg$probe_annot, "Probe annotation table")
annot <- fread(cfg$probe_annot, data.table = FALSE)
if (!cfg$probe_id_col %in% names(annot))
  stop(sprintf("Probe-id column '%s' not in probe annotation.", cfg$probe_id_col))
annot <- annot[match(rownames(expr), annot[[cfg$probe_id_col]]), ]
log_msg("Probe annotation aligned: %d/%d probes matched",
        sum(!is.na(annot[[cfg$probe_id_col]])), nrow(expr))

# ----------------------------------------------------------------------------
# 4. log2 sanity check
# ----------------------------------------------------------------------------
rng <- range(expr, na.rm = TRUE)
is_log2 <- if (!is.na(cfg$already_log2)) cfg$already_log2 else (rng[2] < 100)
log_msg("Value range [%.2f, %.2f] -> treating as %s",
        rng[1], rng[2], ifelse(is_log2, "already log2", "linear intensities"))
if (!is_log2) {
  expr[expr < 1] <- 1           # floor before log to avoid -Inf
  expr <- log2(expr)
  log_msg("Applied log2 transform.")
}

# ----------------------------------------------------------------------------
# 5. Probe-level QC filtering (provider flags; each self-skips if column absent)
# ----------------------------------------------------------------------------
n0 <- nrow(expr)
apply_flag_filter <- function(expr, annot, col, keep_when_true, label) {
  if (is.null(col) || !col %in% names(annot)) {
    log_msg("Probe filter | %s: SKIPPED (column '%s' absent)", label, col); return(expr)
  }
  f <- as_flag(annot[[col]])
  keep <- if (keep_when_true) (f %in% TRUE) else !(f %in% TRUE)
  before <- nrow(expr)
  expr <- expr[keep, , drop = FALSE]
  assign("annot", annot[keep, ], envir = parent.frame())
  log_msg("Probe filter | %s: removed %d (%d remain)", label, before - nrow(expr), nrow(expr))
  expr
}
expr <- apply_flag_filter(expr, annot, cfg$flag_control,   FALSE, "drop control probes")
expr <- apply_flag_filter(expr, annot, cfg$flag_qcok,      TRUE,  "keep probe_QCok")
expr <- apply_flag_filter(expr, annot, cfg$flag_expressed, TRUE,  "keep expressed_blood")
log_msg("Probe QC: %d -> %d probes (%.1f%% retained)",
        n0, nrow(expr), 100 * nrow(expr) / n0)

# Note which retained probes still carry a strong batch association (for the record).
present_assoc <- intersect(cfg$batch_assoc_cols, names(annot))
if (length(present_assoc)) {
  flagged <- Reduce(`|`, lapply(present_assoc, function(c) as_flag(annot[[c]]) %in% TRUE))
  log_msg("Note: %d retained probes flagged with a strong batch association.",
          sum(flagged, na.rm = TRUE))
}

# ----------------------------------------------------------------------------
# 6. Subset to n=48 paired samples, relabel columns -> Sample_Name
# ----------------------------------------------------------------------------
expr_ids   <- link[[cfg$expr_id_col]]
found      <- expr_ids %in% colnames(expr)
if (!all(found))
  log_msg("WARNING: %d paired sampleID(s) absent from expression matrix: %s",
          sum(!found), paste(expr_ids[!found], collapse = ", "))
link       <- link[found, ]
expr_paired <- expr[, link[[cfg$expr_id_col]], drop = FALSE]
colnames(expr_paired) <- link[[cfg$link_meth_name]]   # relabel sampleID -> Sample_Name
log_msg("Paired expression matrix: %d probes x %d samples", nrow(expr_paired), ncol(expr_paired))

pheno_expr <- link
rownames(pheno_expr) <- pheno_expr[[cfg$link_meth_name]]

# ----------------------------------------------------------------------------
# 7. (optional) ComBat batch correction  -- OFF by default
# ----------------------------------------------------------------------------
if (cfg$run_combat) {
  if (!requireNamespace("sva", quietly = TRUE))
    stop("run_combat=TRUE but 'sva' not installed.")
  if (!cfg$combat_batch %in% names(pheno_expr))
    stop("combat_batch column not in pheno: ", cfg$combat_batch)
  b <- factor(pheno_expr[[cfg$combat_batch]])
  if (nlevels(b) < 2 || any(table(b) < 2)) {
    log_msg("ComBat SKIPPED: batch '%s' has <2 usable levels in n=48.", cfg$combat_batch)
  } else {
    log_msg("Running ComBat on batch '%s' (%d levels)...", cfg$combat_batch, nlevels(b))
    expr_paired <- sva::ComBat(dat = expr_paired, batch = b)
  }
} else {
  log_msg("ComBat OFF (expression data ships batch-corrected per DataDictionary).")
}

# ----------------------------------------------------------------------------
# 8. (optional) gene-level collapse
# ----------------------------------------------------------------------------
expr_gene <- NULL
if (cfg$collapse_genes && cfg$gene_symbol_col %in% names(annot)) {
  annot_f <- annot[match(rownames(expr_paired), annot[[cfg$probe_id_col]]), ]
  sym <- annot_f[[cfg$gene_symbol_col]]
  ok  <- !is.na(sym) & sym != "" & sym != "NA"
  # Representative probe per gene: prefer provider's bestexpressed flag, else highest mean.
  if (!is.null(cfg$flag_bestprobe) && cfg$flag_bestprobe %in% names(annot_f)) {
    best <- as_flag(annot_f[[cfg$flag_bestprobe]]) %in% TRUE
    pick <- ok & best
    log_msg("Gene collapse: using provider '%s' (%d probes).", cfg$flag_bestprobe, sum(pick))
  } else {
    mn   <- rowMeans(expr_paired, na.rm = TRUE)
    ord  <- order(sym, -mn)
    pick <- logical(nrow(expr_paired))
    pick[ord[!duplicated(sym[ord])]] <- TRUE
    pick <- pick & ok
    log_msg("Gene collapse: highest-mean probe per gene (%d genes).", sum(pick))
  }
  expr_gene <- expr_paired[pick, , drop = FALSE]
  rownames(expr_gene) <- sym[pick]
  log_msg("Gene-level matrix: %d genes x %d samples", nrow(expr_gene), ncol(expr_gene))
}

# ----------------------------------------------------------------------------
# 9. QC figures + save
# ----------------------------------------------------------------------------
pdf(file.path(cfg$qc_dir, "qc_expr_density.pdf"), width = 7, height = 5)
plot(density(expr_paired[, 1], na.rm = TRUE), main = "Expression density (paired, log2)",
     xlab = "log2 intensity", ylim = c(0, 0.6))
for (j in 2:ncol(expr_paired)) lines(density(expr_paired[, j], na.rm = TRUE), col = "grey60")
dev.off()

# PCA colored by processing batch (visual batch-effect check)
if (nrow(expr_paired) > 2) {
  pc <- prcomp(t(expr_paired), scale. = FALSE)
  bcol <- if ("processing_batch" %in% names(pheno_expr))
            as.integer(factor(pheno_expr$processing_batch)) else 1
  pdf(file.path(cfg$qc_dir, "qc_expr_pca.pdf"), width = 6, height = 6)
  plot(pc$x[, 1], pc$x[, 2], col = bcol, pch = 19,
       xlab = sprintf("PC1 (%.1f%%)", 100 * summary(pc)$importance[2, 1]),
       ylab = sprintf("PC2 (%.1f%%)", 100 * summary(pc)$importance[2, 2]),
       main = "Expression PCA (paired n=48), colored by processing_batch")
  dev.off()
}

annot_paired <- annot[match(rownames(expr_paired), annot[[cfg$probe_id_col]]), ]
saveRDS(expr_paired,  file.path(cfg$out_dir, "expr_probe_paired.rds"))
saveRDS(annot_paired, file.path(cfg$out_dir, "expr_probe_annotation.rds"))
saveRDS(pheno_expr,   file.path(cfg$out_dir, "pheno_expression.rds"))
write.csv(pheno_expr, file.path(cfg$out_dir, "pheno_expression.csv"), row.names = FALSE)
if (!is.null(expr_gene))
  saveRDS(expr_gene, file.path(cfg$out_dir, "expr_gene_paired.rds"))

log_msg("Saved outputs to %s/ and QC figures to %s/", cfg$out_dir, cfg$qc_dir)
writeLines(capture.output(sessionInfo()),
           file.path(cfg$qc_dir, "sessionInfo_02_preprocess_expression.txt"))
log_msg("Done.")
