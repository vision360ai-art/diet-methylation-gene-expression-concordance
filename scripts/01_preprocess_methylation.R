# =============================================================================
# 01_preprocess_methylation.R
#
# Preprocessing of Illumina Infinium MethylationEPIC data (blood, n=98) for the
# LIFE-Adult diet / methylation / expression concordance project.
#
# Pipeline:
#   1. Load & clean the sample sheet
#   2. Read raw IDATs                       -> RGChannelSet
#   3. Sample-level QC (detection p-values, minfi QC, predicted vs reported sex)
#   4. Drop failed samples
#   5. Houseman cell-type deconvolution     -> 6 blood proportions into pheno
#   6. Functional normalization             -> GenomicRatioSet
#   7. Probe-level filtering (detP, SNPs, cross-reactive, multi-hit, [sex chr])
#   8. Extract Beta / M value matrices
#   9. Save processed objects + QC figures
#
# Inputs  : data/raw/SampleSheet_V2.csv, data/raw/idat/*_{Grn,Red}.idat[.gz]
# Outputs : data/processed/*.rds, results/qc/*.pdf|png
#
# Run non-interactively:  Rscript scripts/01_preprocess_methylation.R
# =============================================================================

suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(data.table)
})
# Cell-type deconvolution reference is optional at load time: the Houseman step
# (section 5) self-skips with a warning if this package isn't installed, so the
# rest of preprocessing can still run.
have_flowsorted <- requireNamespace("FlowSorted.Blood.EPIC", quietly = TRUE)
if (have_flowsorted) suppressPackageStartupMessages(library(FlowSorted.Blood.EPIC))

# ----------------------------------------------------------------------------
# 0. Configuration  (edit these; everything below is derived)
# ----------------------------------------------------------------------------
cfg <- list(
  sample_sheet = "data/raw/SampleSheet_V2.csv",
  idat_dir     = "data/raw/idat",     # place *_Grn.idat / *_Red.idat pairs here
  out_dir      = "data/processed",
  qc_dir       = "results/qc",

  # Optional: restrict to methylation+expression PAIRED samples (lower memory, and
  # matches the concordance focus). Uses the expression sample_info for linkage
  # (Synonym == sampleID). Set FALSE to process ALL methylation samples with IDATs.
  restrict_to_paired = TRUE,
  paired_sample_info = "data/raw/Gx_sample_info-21-06-15_v3.txt",

  # --- QC / filtering thresholds ---
  detp_threshold       = 0.01,   # a probe/sample "fails" a well if detP >= this
  sample_detp_max_mean = 0.01,   # drop a sample if its mean detP exceeds this
  probe_detp_max_frac  = 0.01,   # drop a probe failing in > this fraction of samples
  min_idat_bytes       = 1e5,    # per-channel floor; real EPIC IDATs are ~13.6 MB,
                                 # so this only catches empty/truncated/stub files

  # --- switches ---
  # noob_quantile is the default here (not funnorm): funnorm runs noob internally
  # then fits an additional cross-sample regression, giving it the highest peak
  # memory of any step -- it OOMs on this RAM-constrained host. noob+quantile has a
  # strictly lower peak AND is the better-matched normalization for blood diet EWAS,
  # where effects are small/local rather than the large global shifts funnorm targets
  # (Fortin et al. 2014). funnorm remains available as a later sensitivity check.
  normalization        = "noob_quantile",  # "noob_quantile" (default) | "funnorm"
  drop_sex_chromosomes = TRUE,        # blood/diet analysis: usually TRUE
  drop_snp_probes      = TRUE,
  drop_crossreactive   = TRUE,

  seed = 1234
)

set.seed(cfg$seed)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$qc_dir,  recursive = TRUE, showWarnings = FALSE)

# Read IDATs serially. minfi's read.metharray uses BiocParallel internally; the
# Windows default (SnowParam) spawns worker processes that EACH load the ~300 MB
# EPIC manifest, which OOMs on large sample sets. SerialParam reads one at a time.
BiocParallel::register(BiocParallel::SerialParam())

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))

# ----------------------------------------------------------------------------
# 1. Load & clean the sample sheet
# ----------------------------------------------------------------------------
# NOTE: SampleSheet_V2.csv has TWO columns literally named "HLS":
#   - a categorical one ("low"/"high")  -> renamed HLS_group
#   - a numeric composite score          -> renamed HLS_score
# read.csv(check.names = TRUE) would silently rename the 2nd to "HLS.1";
# we handle it explicitly so downstream code is unambiguous.
log_msg("Loading sample sheet: %s", cfg$sample_sheet)
targets <- read.csv(cfg$sample_sheet, check.names = FALSE, stringsAsFactors = FALSE)

hls_cols <- which(names(targets) == "HLS")
if (length(hls_cols) == 2) {
  names(targets)[hls_cols[1]] <- "HLS_group"   # low / high
  names(targets)[hls_cols[2]] <- "HLS_score"   # numeric
} else if ("HLS" %in% names(targets)) {
  names(targets)[names(targets) == "HLS"] <- "HLS_group"
}

# Type coercions
num_cols <- c("Age", "DietScore", "PA_Score", "Smoking_score", "AlcScore",
              "HLS_score", "BMI", "BMI_WAIST_HIP_RATIO", "Weight_in_KG")
for (nc in intersect(num_cols, names(targets))) {
  targets[[nc]] <- suppressWarnings(as.numeric(targets[[nc]]))
}
targets$Sex          <- factor(targets$Gender)                 # "F" / "M"
targets$Sample_Group <- factor(targets$Sample_Group)           # e.g. "4.low"
targets$HLS_group    <- factor(targets$HLS_group, levels = c("low", "high"))

# Basename = absolute path stem minfi appends _Grn.idat / _Red.idat to.
# The sheet's Basename column is "<Sentrix_ID>_<Sentrix_Position>".
targets$Sentrix_ID       <- as.character(targets$Sentrix_ID)
targets$Sentrix_Position <- as.character(targets$Sentrix_Position)
targets$Basename <- file.path(
  normalizePath(cfg$idat_dir, mustWork = FALSE),
  paste0(targets$Sentrix_ID, "_", targets$Sentrix_Position)
)
rownames(targets) <- targets$Sample_Name

log_msg("Sample sheet: %d samples, %d columns", nrow(targets), ncol(targets))
log_msg("  Sample_Group: %s",
        paste(sprintf("%s=%d", levels(targets$Sample_Group),
                      table(targets$Sample_Group)), collapse = ", "))

# --- Optional restriction to paired (methylation+expression) samples ---------
if (isTRUE(cfg$restrict_to_paired) && file.exists(cfg$paired_sample_info)) {
  expr_ids <- data.table::fread(cfg$paired_sample_info)$sampleID
  is_paired <- targets$Synonym %in% expr_ids
  log_msg("Restricting to paired samples (Synonym in expression sampleID): %d/%d kept.",
          sum(is_paired), nrow(targets))
  targets <- targets[is_paired, , drop = FALSE]
}

# --- IDAT presence + integrity filter ---------------------------------------
# The download can be partial (this cohort arrived truncated, e.g. one channel of
# a pair missing at the cut-off point). Rather than fail-fast on the first gap,
# keep only samples whose BOTH channels exist and are non-trivially sized, and log
# every exclusion with a reason for reproducibility.
resolve <- function(stem, chan) {          # prefer plain .idat, fall back to .gz
  p <- paste0(stem, "_", chan, ".idat")
  ifelse(file.exists(p), p, paste0(p, ".gz"))
}
grn <- resolve(targets$Basename, "Grn")
red <- resolve(targets$Basename, "Red")
size_ok <- function(p) { s <- file.size(p); !is.na(s) & s >= cfg$min_idat_bytes }
grn_ok <- file.exists(grn) & size_ok(grn)
red_ok <- file.exists(red) & size_ok(red)
keep_idat <- grn_ok & red_ok

reason <- rep(NA_character_, nrow(targets))
reason[!file.exists(grn) & !file.exists(red)] <- "both IDATs missing"
reason[ file.exists(grn) & !file.exists(red)] <- "Red IDAT missing"
reason[!file.exists(grn) &  file.exists(red)] <- "Grn IDAT missing"
reason[is.na(reason) & !keep_idat]            <- "IDAT below min size (truncated?)"

if (any(!keep_idat)) {
  ex <- data.frame(Sample_Name = targets$Sample_Name[!keep_idat],
                   Basename     = basename(targets$Basename[!keep_idat]),
                   reason       = reason[!keep_idat],
                   stringsAsFactors = FALSE)
  log_msg("Excluding %d/%d sample(s) without a usable IDAT pair:",
          nrow(ex), nrow(targets))
  for (i in seq_len(nrow(ex)))
    log_msg("  - %-8s (%s): %s", ex$Sample_Name[i], ex$Basename[i], ex$reason[i])
  write.csv(ex, file.path(cfg$qc_dir, "excluded_samples_idat.csv"), row.names = FALSE)
}
targets <- targets[keep_idat, , drop = FALSE]
if (!nrow(targets))
  stop("No samples have a complete IDAT pair under '", cfg$idat_dir,
       "'. Check the download or cfg$idat_dir.")
log_msg("Proceeding with %d sample(s) that have complete, non-trivial IDAT pairs.",
        nrow(targets))

# ----------------------------------------------------------------------------
# 2. Read raw IDATs -> RGChannelSet
# ----------------------------------------------------------------------------
log_msg("Reading %d IDAT pairs (this can take a few minutes)...", nrow(targets))
rgSet <- read.metharray.exp(targets = targets, force = TRUE)
pData(rgSet) <- DataFrame(targets)   # attach full phenotype table
sampleNames(rgSet) <- targets$Sample_Name
log_msg("RGChannelSet: %d probes x %d samples | array: %s",
        nrow(rgSet), ncol(rgSet), annotation(rgSet)["array"])

saveRDS(rgSet, file.path(cfg$out_dir, "rgSet_raw.rds"))

# ----------------------------------------------------------------------------
# 3. Sample-level QC
# ----------------------------------------------------------------------------
log_msg("Computing detection p-values...")
detP <- detectionP(rgSet)
sample_mean_detp <- colMeans(detP)

# 3a. Mean detection-p barplot
pdf(file.path(cfg$qc_dir, "qc_mean_detectionP.pdf"), width = 12, height = 5)
barplot(sample_mean_detp, las = 2, cex.names = 0.5,
        ylab = "Mean detection p-value",
        col = ifelse(sample_mean_detp > cfg$sample_detp_max_mean, "red", "grey40"),
        main = "Mean detection p-value per sample")
abline(h = cfg$sample_detp_max_mean, col = "red", lty = 2)
dev.off()

# 3b. minfi's standard QC report (control probes) + getQC (meth vs unmeth medians)
qcReport(rgSet, sampNames = targets$Sample_Name, sampGroups = targets$Sample_Group,
         pdf = file.path(cfg$qc_dir, "qc_report_minfi.pdf"))

mSetRaw <- preprocessRaw(rgSet)
qc <- getQC(mSetRaw)
pdf(file.path(cfg$qc_dir, "qc_meth_vs_unmeth.pdf"), width = 6, height = 6)
plotQC(qc)
dev.off()

# 3c. Predicted vs reported sex (catches sample swaps / mislabels)
gmSet <- mapToGenome(mSetRaw)
pred_sex <- getSex(gmSet)
targets$predictedSex <- pred_sex$predictedSex
sex_mismatch <- toupper(substr(targets$predictedSex, 1, 1)) !=
                toupper(substr(as.character(targets$Sex), 1, 1))
# Plot manually rather than minfi::plotSex(): that helper calls colData() on the
# getSex DFrame, which fails under this minfi/S4Vectors version. Same diagnostic.
pdf(file.path(cfg$qc_dir, "qc_predicted_sex.pdf"), width = 6, height = 6)
plot(pred_sex$xMed, pred_sex$yMed,
     col = ifelse(pred_sex$predictedSex == "M", "steelblue", "firebrick"),
     pch = 16, xlab = "chrX median total intensity (log2)",
     ylab = "chrY median total intensity (log2)", main = "Predicted sex (getSex)")
text(pred_sex$xMed, pred_sex$yMed, labels = targets$Sample_Name, cex = 0.5, pos = 3)
legend("bottomleft", c("predicted M", "predicted F"),
       col = c("steelblue", "firebrick"), pch = 16, bty = "n")
dev.off()
if (any(sex_mismatch, na.rm = TRUE)) {
  log_msg("WARNING: %d sample(s) with predicted != reported sex: %s",
          sum(sex_mismatch, na.rm = TRUE),
          paste(targets$Sample_Name[which(sex_mismatch)], collapse = ", "))
}

# Free the large QC intermediates (raw MethylSet + genome-mapped set) before the
# memory-heavy normalization -- they are not needed past the sex check. On a
# RAM-constrained host this is the difference between funnorm fitting or OOMing.
rm(mSetRaw, gmSet, qc, pred_sex); invisible(gc())

# ----------------------------------------------------------------------------
# 4. Drop failed samples
# ----------------------------------------------------------------------------
keep_samples <- sample_mean_detp <= cfg$sample_detp_max_mean
log_msg("Sample QC: keeping %d / %d (dropped: %s)",
        sum(keep_samples), length(keep_samples),
        ifelse(all(keep_samples), "none",
               paste(names(keep_samples)[!keep_samples], collapse = ", ")))

rgSet   <- rgSet[, keep_samples]
detP    <- detP[, keep_samples]
targets <- targets[keep_samples, ]

# ----------------------------------------------------------------------------
# 5. Houseman cell-type deconvolution (blood composition)
# ----------------------------------------------------------------------------
# Blood cell composition is the dominant confounder for whole-blood methylation.
# estimateCellCounts2() (FlowSorted.Blood.EPIC) reconstructs each sample's
# proportion of 6 leukocyte subsets by co-normalizing against a sorted-blood
# reference. It MUST run on the raw RGChannelSet (before normalization/filtering)
# because it needs the reference CpGs at native intensities. We use the
# IDOL-optimized CpG library, the current best practice for adult EPIC blood.
#
# The 6 proportions are written as columns into pData(rgSet) so they propagate
# through preprocessFunnorm into grSet's pheno and become model covariates
# downstream (no re-read, no join needed).
#
# Memory note: estimateCellCounts2 co-normalizes with a ~37-sample reference
# (~73 samples x 1.05M probes, noob) and OOMs when run inline here on a RAM-
# constrained host (it fires while rgSet + detP are both resident). So we PREFER a
# cached table from scripts/01a_cell_counts.R -- an isolated process that holds only
# rgSet + the reference. If the cache is absent we fall back to the inline call and
# write the cache for next time.
cell_cache <- file.path(cfg$out_dir, "cell_counts.csv")
cell_counts <- NULL

if (file.exists(cell_cache)) {
  log_msg("Loading cached cell-type proportions: %s", cell_cache)
  cc <- read.csv(cell_cache, row.names = 1, check.names = FALSE)
  missing <- setdiff(sampleNames(rgSet), rownames(cc))
  if (length(missing)) {
    log_msg("  WARNING: cache is missing %d sample(s): %s -- ignoring cache, recomputing.",
            length(missing), paste(missing, collapse = ", "))
  } else {
    cell_counts <- cc[sampleNames(rgSet), , drop = FALSE]   # enforce order
  }
}

if (is.null(cell_counts) && have_flowsorted) {
  log_msg("Estimating blood cell-type proportions inline (Houseman / IDOL)...")
  log_msg("  NOTE: high peak memory -- if this OOMs, run scripts/01a_cell_counts.R")
  log_msg("        (isolated process) to produce %s, then re-run 01.", cell_cache)
  cell_est <- estimateCellCounts2(
    rgSet,
    compositeCellType   = "Blood",
    processMethod       = "preprocessNoob",
    probeSelect         = "IDOL",
    cellTypes           = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
    referencePlatform   = "IlluminaHumanMethylationEPIC",
    IDOLOptimizedCpGs   = IDOLOptimizedCpGsBloodEPIC,
    returnAll           = FALSE
  )
  cell_counts <- as.data.frame(cell_est$counts)     # samples x 6, rows = sample names
  cell_counts <- cell_counts[sampleNames(rgSet), , drop = FALSE]  # enforce order
  write.csv(cell_counts, cell_cache, row.names = TRUE)            # cache for re-runs
  log_msg("  Cached cell proportions -> %s", cell_cache)
}

if (!is.null(cell_counts)) {
  for (ct in colnames(cell_counts)) {
    pData(rgSet)[[ct]] <- cell_counts[[ct]]
    targets[[ct]]      <- cell_counts[[ct]]
  }
  log_msg("  Cell proportions (mean): %s",
          paste(sprintf("%s=%.3f", colnames(cell_counts),
                        colMeans(cell_counts)), collapse = ", "))

  # Cell-composition QC plot
  pdf(file.path(cfg$qc_dir, "qc_cell_composition.pdf"), width = 10, height = 5)
  boxplot(cell_counts, ylab = "Estimated proportion", las = 2,
          main = "Houseman blood cell-type proportions", col = "grey80")
  dev.off()
} else {
  log_msg("SKIPPING Houseman cell-type step: no cache and FlowSorted.Blood.EPIC not usable.")
  log_msg("  -> pheno will lack CD8T/CD4T/NK/Bcell/Mono/Neu columns. Install the")
  log_msg("     package (scripts/install_flowsorted.R) and run scripts/01a_cell_counts.R,")
  log_msg("     then re-run 01; 03's design self-adjusts to whichever covariates exist.")
}

# ----------------------------------------------------------------------------
# 6. Normalization -> GenomicRatioSet
# ----------------------------------------------------------------------------
log_msg("Normalizing (method = %s)...", cfg$normalization)
if (cfg$normalization == "funnorm") {
  # Functional normalization: recommended when global methylation differences
  # are expected between groups (here: lifestyle-extreme case/control design).
  # NOTE: highest peak memory of any step -- can OOM on RAM-constrained hosts.
  grSet <- preprocessFunnorm(rgSet)
} else if (cfg$normalization == "noob_quantile") {
  # Lower-peak path. Free the raw RGChannelSet between the two stages: nothing
  # downstream needs it once normalization starts (Houseman + detP already ran,
  # probe filtering works off grSet + the cached detP), and dropping it here is
  # what keeps quantile normalization inside the memory envelope.
  mSet  <- preprocessNoob(rgSet)
  rm(rgSet); invisible(gc())
  grSet <- preprocessQuantile(mSet)   # returns a GenomicRatioSet
  rm(mSet); invisible(gc())
} else {
  stop("Unknown cfg$normalization: ", cfg$normalization)
}
log_msg("Normalized GenomicRatioSet: %d probes x %d samples",
        nrow(grSet), ncol(grSet))

# ----------------------------------------------------------------------------
# 7. Probe-level filtering
# ----------------------------------------------------------------------------
# Align detP to the normalized object (same probe order & samples).
detP <- detP[match(rownames(grSet), rownames(detP)), colnames(grSet), drop = FALSE]

n_start <- nrow(grSet)

# 7a. Probes failing detection in too many samples
frac_failed <- rowMeans(detP >= cfg$detp_threshold)
keep_probes <- frac_failed <= cfg$probe_detp_max_frac
grSet <- grSet[keep_probes, ]
log_msg("Probe filter | detection-p: removed %d (%.1f%%)",
        n_start - nrow(grSet), 100 * (n_start - nrow(grSet)) / n_start)

# 7b. SNP-affected probes (SBE / CpG SNPs)
if (cfg$drop_snp_probes) {
  before <- nrow(grSet)
  grSet  <- dropLociWithSnps(grSet)
  log_msg("Probe filter | SNP loci: removed %d", before - nrow(grSet))
}

# 7c. Sex-chromosome probes
if (cfg$drop_sex_chromosomes) {
  before <- nrow(grSet)
  ann    <- getAnnotation(grSet)
  autosome <- !(ann$chr %in% c("chrX", "chrY"))
  grSet  <- grSet[autosome, ]
  log_msg("Probe filter | sex chromosomes: removed %d", before - nrow(grSet))
}

# 7d. Cross-reactive / multi-hit probes (Pidsley et al. 2016, EPIC).
#     Provide the annotation CSV at data/raw/annotation/Pidsley2016_crossreactive.csv
#     (one probe ID per line, or a column named "TargetID"/"IlmnID"/"ProbeID").
if (cfg$drop_crossreactive) {
  xr_path <- "data/raw/annotation/Pidsley2016_crossreactive.csv"
  if (file.exists(xr_path)) {
    xr <- fread(xr_path, header = TRUE)
    id_col <- intersect(c("TargetID", "IlmnID", "ProbeID", "V1"), names(xr))[1]
    xr_ids <- unique(as.character(xr[[id_col]]))
    before <- nrow(grSet)
    grSet  <- grSet[!(rownames(grSet) %in% xr_ids), ]
    log_msg("Probe filter | cross-reactive: removed %d", before - nrow(grSet))
  } else {
    log_msg("Probe filter | cross-reactive: SKIPPED (list not found at %s)",
            xr_path)
  }
}

log_msg("Final probe set: %d (%.1f%% of %d retained)",
        nrow(grSet), 100 * nrow(grSet) / n_start, n_start)

# ----------------------------------------------------------------------------
# 8. Extract Beta / M value matrices
# ----------------------------------------------------------------------------
beta <- getBeta(grSet)
# Guard against +/-Inf M-values from beta at exactly 0 or 1, then derive M.
beta[beta < 1e-6]     <- 1e-6
beta[beta > 1 - 1e-6] <- 1 - 1e-6
mval <- log2(beta / (1 - beta))

log_msg("Beta matrix: %d x %d | M matrix: %d x %d",
        nrow(beta), ncol(beta), nrow(mval), ncol(mval))

# 8a. Beta density (probe-type / normalization sanity check)
pdf(file.path(cfg$qc_dir, "qc_beta_density_normalized.pdf"), width = 7, height = 5)
densityPlot(beta, sampGroups = targets$Sample_Group,
            main = "Beta density (normalized, filtered)")
dev.off()

# ----------------------------------------------------------------------------
# 9. Save processed objects
# ----------------------------------------------------------------------------
pheno <- as.data.frame(pData(grSet))   # includes the 6 cell-type proportions

# Robust writer: on this host, antivirus real-time scanning intermittently locks a
# freshly-written large .rds and makes saveRDS fail with "error writing to
# connection". Retry a few times (removing the partial file first) so a transient
# lock doesn't discard a completed normalization.
robust_saveRDS <- function(obj, path, tries = 4, wait = 5) {
  for (i in seq_len(tries)) {
    ok <- tryCatch({ saveRDS(obj, path); TRUE },
                   error = function(e) { log_msg("  saveRDS(%s) attempt %d/%d failed: %s",
                                                 basename(path), i, tries, conditionMessage(e)); FALSE })
    if (ok) return(invisible(TRUE))
    if (file.exists(path)) unlink(path)      # drop the truncated/0-byte file
    Sys.sleep(wait); invisible(gc())
  }
  stop("Could not write ", path, " after ", tries, " attempts.")
}

robust_saveRDS(grSet, file.path(cfg$out_dir, "grSet_normalized_filtered.rds"))
robust_saveRDS(beta,  file.path(cfg$out_dir, "beta_values.rds"))
robust_saveRDS(mval,  file.path(cfg$out_dir, "m_values.rds"))
robust_saveRDS(pheno, file.path(cfg$out_dir, "pheno_methylation.rds"))
write.csv(pheno, file.path(cfg$out_dir, "pheno_methylation.csv"), row.names = FALSE)

log_msg("Saved outputs to %s/ and QC figures to %s/", cfg$out_dir, cfg$qc_dir)

# Reproducibility record
writeLines(capture.output(sessionInfo()),
           file.path(cfg$qc_dir, "sessionInfo_01_preprocess_methylation.txt"))
log_msg("Done.")
