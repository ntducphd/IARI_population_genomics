#!/usr/bin/env Rscript
# 17_phase_haplo_scans.R -- [stage 17] haplotype-based selection scans on Set 1:
# Beagle phasing -> rehh iHS (panel-wide) and XP-EHH (between the two largest admixture clusters).
#
# Motivation: the modern, haplotype-length selection
# layer. Set 1 is whole-genome-resequenced (502,675 QC SNPs), so phased-haplotype statistics are
# feasible; Set 2 (50K chip) is excluded (marker density far below EHH requirements).
#
# Honest caveat carried into the manuscript: EHH-family statistics were developed for outcrossing
# populations; near-complete selfing (F ~ 0.93, stage 16) lengthens haplotypes genome-wide and
# reduces effective recombination, so absolute iHS/XP-EHH values are not calibrated against the
# outcrossing null -- we therefore interpret only the RELATIVE ranking of regions (empirical
# outliers, |standardised score| in the top 0.5% / >= 3), cross-checked against the pcadapt scan
# (stage 18), and state this limitation explicitly. Phasing itself is nearly deterministic here:
# with ~93% of the genome homozygous, few sites are genuinely ambiguous.
#
# Pipeline:
#   [1] PLINK Set1_qc bfile -> per run: use existing Set1_qc.vcf (stage 15 recode)
#   [2] Beagle 5.5 (analysis/tools/beagle.jar, Java 17): phase, impute=false
#   [3] rehh per chromosome: data2haplohh -> scan_hh -> ihh2ihs (iHS), ies2xpehh (XP-EHH,
#       cluster C-largest vs C-second from the stage-03 admixture argmax assignment)
#   [4] outputs: per-SNP scores + empirical-outlier flags + a candidate-region summary
#      (adjacent outlier SNPs within 200 kb merged)
#
# Outputs (analysis/results/tables/):
#   ihs_set1.csv            -- chrom, pos, ihs, p_bilateral, outlier_top (|iHS|>=3)
#   xpehh_set1.csv          -- chrom, pos, xpehh, p_bilateral, outlier_top (|XP-EHH|>=3)
#   selection_regions_set1.csv -- merged outlier regions x {iHS, XPEHH, pcadapt overlap}
#   haplo_scan_summary.csv  -- n SNPs scanned, n outliers per statistic, cluster sizes used
# NOTE: run via PowerShell (house SVD rule); long-running -- launched in background.

suppressMessages({ library(rehh); library(vcfR) })

args      <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPTS   <- dirname(normalizePath(this_file))
ROOT      <- normalizePath(file.path(SCRIPTS, "..", ".."))
INTERIM   <- file.path(ROOT, "analysis", "data", "interim")
TAB       <- file.path(ROOT, "analysis", "results", "tables")
TOOLS     <- file.path(ROOT, "analysis", "tools")

phased_vcf <- file.path(INTERIM, "Set1_phased.vcf.gz")

# [2] Beagle phasing (skipped if output already present)
if (!file.exists(phased_vcf)) {
  cmd <- sprintf('java -Xmx6g -jar "%s" gt="%s" out="%s" impute=false nthreads=4',
                 file.path(TOOLS, "beagle.jar"),
                 file.path(INTERIM, "Set1_qc.vcf"),
                 file.path(INTERIM, "Set1_phased"))
  cat("[17] phasing with Beagle...\n")
  status <- system(cmd)
  if (status != 0) stop("Beagle phasing failed")
}

# cluster assignment (argmax Q) for XP-EHH populations
q      <- read.csv(file.path(TAB, "admixture_set1_Q.csv"))
qcols  <- grep("^Q", names(q), value = TRUE)
assign <- paste0("C", max.col(q[qcols]) - 1L)
names(assign) <- q$sample
sizes  <- sort(table(assign), decreasing = TRUE)
popA   <- names(sizes)[1]; popB <- names(sizes)[2]
cat(sprintf("[17] XP-EHH populations: %s (n=%d) vs %s (n=%d)\n",
            popA, sizes[1], popB, sizes[2]))

samplesA <- names(assign)[assign == popA]
samplesB <- names(assign)[assign == popB]

chroms <- as.character(1:12)
ihs_scan <- NULL; xp_scan <- NULL

for (ch in chroms) {
  hhA <- tryCatch(
    data2haplohh(hap_file = phased_vcf, chr.name = ch, polarize_vcf = FALSE,
                 vcf_reader = "vcfR", verbose = FALSE),
    error = function(e) NULL)
  if (is.null(hhA)) { cat(sprintf("[17] chr %s: no data, skipped\n", ch)); next }

  scanAll <- scan_hh(hhA, polarized = FALSE, discard_integration_at_border = TRUE)
  ihs_scan <- rbind(ihs_scan, scanAll)

  hh_a <- subset(hhA, select.hap = which(sub("_[12]$", "", hap.names(hhA)) %in% samplesA))
  hh_b <- subset(hhA, select.hap = which(sub("_[12]$", "", hap.names(hhA)) %in% samplesB))
  sa <- scan_hh(hh_a, polarized = FALSE, discard_integration_at_border = TRUE)
  sb <- scan_hh(hh_b, polarized = FALSE, discard_integration_at_border = TRUE)
  xp <- ies2xpehh(sa, sb, popname1 = popA, popname2 = popB, verbose = FALSE)
  xp_scan <- rbind(xp_scan, xp)
  cat(sprintf("[17] chr %s: %d SNPs scanned\n", ch, nrow(scanAll)))
}

# iHS (unpolarized -> uses |iHS|; frequency-bin standardisation by ihh2ihs)
ihs <- ihh2ihs(ihs_scan, freqbin = 0.05, verbose = FALSE)
ihs_df <- data.frame(chrom = ihs$ihs$CHR, pos = ihs$ihs$POSITION,
                     ihs = ihs$ihs$IHS, logp = ihs$ihs$LOGPVALUE)
ihs_df$outlier <- !is.na(ihs_df$ihs) & abs(ihs_df$ihs) >= 3
write.csv(ihs_df, file.path(TAB, "ihs_set1.csv"), row.names = FALSE)

xp_df <- data.frame(chrom = xp_scan$CHR, pos = xp_scan$POSITION,
                    xpehh = xp_scan[[grep("XPEHH", names(xp_scan), value = TRUE)[1]]],
                    logp = xp_scan[[grep("LOGPVALUE", names(xp_scan), value = TRUE)[1]]])
xp_df$outlier <- !is.na(xp_df$xpehh) & abs(xp_df$xpehh) >= 3
write.csv(xp_df, file.path(TAB, "xpehh_set1.csv"), row.names = FALSE)

# merged candidate regions (outlier SNPs within 200 kb merged), + pcadapt overlap if available
merge_regions <- function(df, score_name) {
  d <- df[df$outlier, c("chrom", "pos")]
  if (nrow(d) == 0) return(data.frame())
  d <- d[order(d$chrom, d$pos), ]
  gap <- 2e5
  d$new <- c(TRUE, diff(d$pos) > gap | d$chrom[-1] != d$chrom[-nrow(d)])
  d$region <- cumsum(d$new)
  agg <- aggregate(pos ~ chrom + region, d, function(p) c(min(p), max(p), length(p)))
  data.frame(chrom = agg$chrom, start = agg$pos[, 1], end = agg$pos[, 2],
             n_outlier_snps = agg$pos[, 3], statistic = score_name)
}
regions <- rbind(merge_regions(ihs_df, "iHS"), merge_regions(xp_df, "XP-EHH"))
pc_path <- file.path(TAB, "pcadapt_outliers_set1.csv")
if (file.exists(pc_path) && nrow(regions) > 0) {
  pc <- read.csv(pc_path); pc <- pc[pc$outlier, ]
  regions$pcadapt_overlap <- mapply(function(chn, s, e)
    any(pc$chrom == chn & pc$pos >= s - 1e5 & pc$pos <= e + 1e5),
    regions$chrom, regions$start, regions$end)
}
write.csv(regions, file.path(TAB, "selection_regions_set1.csv"), row.names = FALSE)

write.csv(data.frame(
  n_snps_ihs = sum(!is.na(ihs_df$ihs)), n_outlier_ihs = sum(ihs_df$outlier),
  n_snps_xpehh = sum(!is.na(xp_df$xpehh)), n_outlier_xpehh = sum(xp_df$outlier),
  popA = popA, nA = as.integer(sizes[1]), popB = popB, nB = as.integer(sizes[2]),
  n_regions = nrow(regions)),
  file.path(TAB, "haplo_scan_summary.csv"), row.names = FALSE)
cat(sprintf("[17] done: %d iHS outliers, %d XP-EHH outliers, %d merged regions\n",
            sum(ihs_df$outlier), sum(xp_df$outlier), nrow(regions)))
