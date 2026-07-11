#!/usr/bin/env Rscript
# =============================================================================
# 03b_annotate_hits.R  --  annotate the EWAS discovery hits (nominal p<1e-3) with
# EPIC gene/region info and tally candidate-gene vs. novel.
#
# Reads the tracked trimmed hits table from 03a_ewas.R, maps each CpG to gene(s)
# and region(s) via IlluminaHumanMethylationEPICanno.ilm10b4.hg19, and checks each
# gene against the Klemp candidate_genes.txt (exact UCSC symbol match, same rule
# 03_concordance.R uses to flag candidate pairs).
#
# Input  : results/tables/ewas_diet_significant.csv   (996 CpGs + stats)
#          data/raw/annotation/candidate_genes.txt
# Output : results/tables/ewas_hits_annotated.csv      (per-CpG annotation + flags)
# =============================================================================
user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib) && dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(data.table)
})

hits <- fread("results/tables/ewas_diet_significant.csv")
cat(sprintf("Loaded %d EWAS hits (nominal p<1e-3).\n", nrow(hits)))

ann  <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
cols <- c("chr","pos","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_Island")
a    <- as.data.frame(ann[hits$CpG, cols, drop = FALSE])
hits <- cbind(hits, a)

# per-CpG unique gene list (UCSC_RefGene_Name is ";"-separated, often repeated/transcript)
gene_list <- lapply(hits$UCSC_RefGene_Name, function(x)
  if (is.na(x) || x == "") character(0) else unique(strsplit(x, ";", fixed = TRUE)[[1]]))
hits$n_genes <- lengths(gene_list)
hits$genes   <- vapply(gene_list, paste, "", collapse = "|")

# candidate list (same parse as 03: strip full-line + inline # comments)
raw  <- readLines("data/raw/annotation/candidate_genes.txt")
cand <- unique(trimws(sub("#.*$", "", raw))); cand <- cand[cand != ""]

hits$candidate_genes_hit <- vapply(gene_list, function(g)
  paste(intersect(g, cand), collapse = "|"), "")
hits$on_candidate <- vapply(gene_list, function(g) any(g %in% cand), logical(1))
hits$intergenic   <- hits$n_genes == 0

out_cols <- c("CpG","logFC","P.Value","adj.P.Val","chr","pos","genes",
              "UCSC_RefGene_Group","Relation_to_Island","n_genes",
              "on_candidate","candidate_genes_hit")
fwrite(hits[, ..out_cols], "results/tables/ewas_hits_annotated.csv")

# ---- tallies --------------------------------------------------------------
n_total   <- nrow(hits)
n_inter   <- sum(hits$intergenic)
n_genic   <- n_total - n_inter
n_cpg_cand<- sum(hits$on_candidate)

all_genes  <- sort(unique(unlist(gene_list)))
cand_hit   <- sort(intersect(all_genes, cand))
novel_gene <- setdiff(all_genes, cand)

cat("\n================= CpG-LEVEL TALLY (n=", n_total, ") =================\n", sep="")
cat(sprintf("  intergenic (no gene)          : %4d (%.1f%%)\n", n_inter, 100*n_inter/n_total))
cat(sprintf("  gene-annotated                : %4d (%.1f%%)\n", n_genic, 100*n_genic/n_total))
cat(sprintf("    -> on a candidate gene      : %4d (%.1f%%)\n", n_cpg_cand, 100*n_cpg_cand/n_total))
cat(sprintf("    -> novel (genic, non-cand)  : %4d (%.1f%%)\n", n_genic-n_cpg_cand,
            100*(n_genic-n_cpg_cand)/n_total))

cat("\n================= GENE-LEVEL TALLY =================\n")
cat(sprintf("  unique genes represented      : %d\n", length(all_genes)))
cat(sprintf("  candidate genes hit           : %d of %d on the list\n", length(cand_hit), length(cand)))
cat(sprintf("  novel genes (not on list)     : %d\n", length(novel_gene)))

if (length(cand_hit)) {
  cat("\nCandidate genes with >=1 hit CpG:\n")
  for (g in cand_hit) {
    idx  <- vapply(gene_list, function(gg) g %in% gg, logical(1))
    best <- hits[idx][which.min(P.Value)]
    cat(sprintf("  %-10s  %d CpG(s), best %s p=%.2e (%s)\n", g, sum(idx),
                best$CpG, best$P.Value, best$UCSC_RefGene_Group))
  }
} else cat("\nNo candidate genes among the 996 hits.\n")

cat("\nTop 5 hits overall (gene | region | candidate?):\n")
top <- hits[order(P.Value)][1:5]
for (i in 1:nrow(top))
  cat(sprintf("  %-12s %-10s %-14s cand=%s  p=%.2e\n", top$CpG[i],
              ifelse(top$genes[i]=="", "<intergenic>", top$genes[i]),
              top$UCSC_RefGene_Group[i], top$on_candidate[i], top$P.Value[i]))

cat("\nWrote results/tables/ewas_hits_annotated.csv\nDONE\n")
