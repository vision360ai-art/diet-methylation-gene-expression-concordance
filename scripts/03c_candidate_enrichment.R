#!/usr/bin/env Rscript
# =============================================================================
# 03c_candidate_enrichment.R  --  is the Klemp candidate-gene set over-represented
# among the EWAS discovery hits, beyond chance?
#
# Hypergeometric test at two levels, against the universe of ALL tested CpGs:
#   CpG-level  : are candidate-gene CpGs enriched among the 996 hit CpGs?
#   gene-level : are candidate genes enriched among the unique hit genes? (cleaner --
#                CpGs within a gene are correlated, so CpG-level pseudo-replicates.)
# Candidate matching is exact UCSC symbol (same rule 03_concordance.R uses).
#
# Interpretation caveats recorded with the result: at p<1e-3 over ~807k tests ~807
# hits are expected under the null (n=36 is underpowered), so the discovery pool is
# mostly chance; and this is a single un-corrected test on small counts.
#
# Inputs  : results/tables/ewas_diet_significant.csv   (hit CpGs, from 03a)
#           results/tables/ewas_diet_methylation.csv.gz (full tested set; or falls
#             back to data/processed/m_values.rds rownames) -- both regenerable by 03a
#           data/raw/annotation/candidate_genes.txt
# Output  : results/tables/candidate_enrichment.txt     (the numbers, tracked)
# =============================================================================
user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib) && dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(data.table)
})

cfg <- list(
  hits    = "results/tables/ewas_diet_significant.csv",
  full_gz = "results/tables/ewas_diet_methylation.csv.gz",
  mvals   = "data/processed/m_values.rds",
  cand    = "data/raw/annotation/candidate_genes.txt",
  out     = "results/tables/candidate_enrichment.txt"
)

# report buffer: everything cat'd is also written to cfg$out
REP <- character(0)
say <- function(...) { line <- sprintf(...); cat(line, "\n"); REP[[length(REP)+1]] <<- line }

# --- tested universe of CpGs ----------------------------------------------
if (!file.exists(cfg$hits)) stop("Missing ", cfg$hits, " -- run 03a_ewas.R first.")
hits <- fread(cfg$hits, select = "CpG")$CpG
if (file.exists(cfg$full_gz)) {
  allcpg <- fread(cfg$full_gz, select = "CpG")$CpG
} else if (file.exists(cfg$mvals)) {
  allcpg <- rownames(readRDS(cfg$mvals))
} else stop("Need the full tested CpG set: run 03a_ewas.R (writes the .gz) or keep m_values.rds.")

N <- length(allcpg); n <- length(hits)
raw  <- readLines(cfg$cand)
cand <- unique(trimws(sub("#.*$", "", raw))); cand <- cand[cand != ""]
say("Universe (tested CpGs) N=%d ; discovery hits n=%d ; candidate genes=%d", N, n, length(cand))

# --- map CpGs -> gene sets -------------------------------------------------
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
gene_of <- function(cpgs)
  lapply(strsplit(ann[cpgs, "UCSC_RefGene_Name"], ";", fixed = TRUE),
         function(x) if (length(x) == 0) character(0) else unique(x))
g_all <- gene_of(allcpg); g_hit <- gene_of(hits)

# --- CpG-level hypergeometric ---------------------------------------------
K <- sum(vapply(g_all, function(x) any(x %in% cand), logical(1)))   # candidate CpGs in universe
k <- sum(vapply(g_hit, function(x) any(x %in% cand), logical(1)))   # candidate CpGs among hits
exp_cpg <- n * K / N
say("")
say("===== CpG-LEVEL hypergeometric =====")
say("  candidate CpGs in universe K = %d (%.3f%% of %d)", K, 100*K/N, N)
say("  observed among hits k = %d ; expected = %.2f ; fold = %.2fx", k, exp_cpg, k/exp_cpg)
say("  P(X>=%d) over-representation  = %.3g", k, phyper(k-1, K, N-K, n, lower.tail = FALSE))
say("  P(X<=%d) under-representation = %.3g", k, phyper(k,   K, N-K, n, lower.tail = TRUE))

# --- gene-level hypergeometric --------------------------------------------
univ_genes <- unique(unlist(g_all)); hit_genes <- unique(unlist(g_hit))
G <- length(univ_genes)
C <- length(intersect(cand, univ_genes))            # candidates actually on the array
nn <- length(hit_genes); kk <- length(intersect(cand, hit_genes))
exp_gene <- nn * C / G
say("")
say("===== GENE-LEVEL hypergeometric =====")
say("  gene universe G = %d ; candidates present C = %d of %d (exact symbol)", G, C, length(cand))
miss <- setdiff(cand, univ_genes)
say("  candidates absent from array universe: %s", if (length(miss)) paste(miss, collapse=", ") else "(none)")
say("  hit genes = %d ; candidate hit genes kk = %d ; expected = %.2f ; fold = %.2fx",
    nn, kk, exp_gene, kk/exp_gene)
say("  P(X>=%d) over-representation  = %.3g", kk, phyper(kk-1, C, G-C, nn, lower.tail = FALSE))
say("  P(X<=%d) under-representation = %.3g", kk, phyper(kk,   C, G-C, nn, lower.tail = TRUE))

say("")
say("CAVEATS: n=36 underpowered -> ~%d of %d hits expected under the null at p<1e-3;",
    round(N*1e-3), n)
say("  small counts (kk=%d, k=%d); single un-corrected test -> treat as weak/borderline.", kk, k)

dir.create(dirname(cfg$out), recursive = TRUE, showWarnings = FALSE)
writeLines(REP, cfg$out)
cat("\nWrote", cfg$out, "\nDONE\n")
