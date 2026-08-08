#!/usr/bin/env Rscript
# 18_selection_outliers.R -- [stage 18] Fst-outlier selection scan with pcadapt, both panels.
#
# Motivation: first selection layer of the paper.
# pcadapt detects loci whose allele frequencies are atypically structured along the panel's own
# principal components (Luu et al. 2017) -- no phasing, no external reference, chip-safe, and
# robust to the admixed, selfing composition of these panels (it conditions on the realised
# structure rather than on discrete populations).
#
# OutFLANK was planned as the second outlier method but its hard dependency (Bioconductor
# `qvalue`) is uninstallable on this machine (Bioconductor repository unreachable from the
# analysis environment, 2026-08-08) -- documented fallback per house convention: pcadapt is the
# reported genome-scan; the haplotype-based scans (stage 17) provide the independent second line.
# Multiple-testing control uses Benjamini-Hochberg FDR (stats::p.adjust), q < 0.05.
#
# K choice: the number of PCs retained follows each panel's structure analysis (elbow K from
# stage 03 minus 1 is the natural upper bound for axes of structure; we scan K = 2..10 scree and
# fix K at the elbow of pcadapt's own singular-value scree, reported in the output for audit).
#
# Outputs (analysis/results/tables/):
#   pcadapt_outliers_{set}.csv  -- chrom, pos, snp, pvalue, qvalue_bh, outlier (q<0.05)
#   pcadapt_summary.csv         -- per panel: K_used, n_snps, n_outliers, prop_outliers
# NOTE (run environment, 2026-08-08): must be run via PowerShell on Windows -- Rscript under
# Git Bash/MSYS2 segfaults on SVD-family calls (documented machine constraint; run R stages via PowerShell).

suppressMessages({ library(pcadapt) })

args      <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPTS   <- dirname(normalizePath(this_file))
ROOT      <- normalizePath(file.path(SCRIPTS, "..", ".."))
INTERIM   <- file.path(ROOT, "analysis", "data", "interim")
TAB       <- file.path(ROOT, "analysis", "results", "tables")

summary_rows <- list()

for (panel in c("Set1", "Set2")) {
  bed <- file.path(INTERIM, paste0(panel, "_qc.bed"))
  bim <- read.table(file.path(INTERIM, paste0(panel, "_qc.bim")),
                    col.names = c("chrom", "snp", "cm", "pos", "a1", "a2"))
  dat <- read.pcadapt(bed, type = "bed")

  # scree scan to choose K: largest K before the singular values flatten (Cattell rule as
  # implemented by inspecting proportion of explained variance drops)
  scree <- pcadapt(dat, K = 10, method = "mahalanobis")
  pv    <- scree$singular.values^2
  drops <- -diff(pv) / pv[-length(pv)]
  K     <- max(2L, which(drops < 0.10)[1])          # first K where the next axis adds <10%
  if (is.na(K)) K <- 5L

  fit <- pcadapt(dat, K = K, method = "mahalanobis")
  q   <- p.adjust(fit$pvalues, method = "BH")
  out <- data.frame(chrom = bim$chrom, pos = bim$pos, snp = bim$snp,
                    pvalue = fit$pvalues, qvalue_bh = q,
                    outlier = !is.na(q) & q < 0.05)
  write.csv(out, file.path(TAB, paste0("pcadapt_outliers_", tolower(panel), ".csv")),
            row.names = FALSE)

  n_out <- sum(out$outlier, na.rm = TRUE)
  summary_rows[[panel]] <- data.frame(panel = panel, K_used = K, n_snps = nrow(out),
                                      n_outliers = n_out,
                                      prop_outliers = n_out / nrow(out))
  cat(sprintf("[18] %s: K = %d, %d/%d SNPs outliers at BH q < 0.05 (%.2f%%)\n",
              panel, K, n_out, nrow(out), 100 * n_out / nrow(out)))
}

write.csv(do.call(rbind, summary_rows), file.path(TAB, "pcadapt_summary.csv"),
          row.names = FALSE)
cat("[18] done -> pcadapt_summary.csv\n")
