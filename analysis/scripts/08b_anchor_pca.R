# 08b_anchor_pca.R — project Set1 into the global 3K-RGP reference PCA space (3024 accessions:
# 2874 3K-RGP reference + our 150 Set1, already merged and LD-pruned at pruned_v2.1). Loads ONLY
# SNPRelate + gdsfmt (the tested-safe pair; ape is NOT loaded here, per the established rule).
#
# The full pruned_v2.1 (3024 x 1.01M SNPs) was tried directly and abandoned: BED->GDS conversion
# alone did not finish within ~19 minutes (killed). Thinned to a random ~10% subsample via PLINK
# first (same, already-established pattern as 07_ld_decay.py's Set1 thinning) -- a PCA structure
# plot does not need every one of 1.01M already-pruned markers, and this keeps the whole anchor
# step to a few minutes.
#
# Set2 (50K array) is NOT anchored here: no equivalent pre-merged 3K-RGP+Set2 genotype file was
# found on disk (only a common-SNP-position list, 3K-HDRA-snp-comm-miss5pc.txt, without the actual
# 3K-RGP genotypes at those positions) -- building that anchor would require additional data
# preparation out of scope for this pass. Documented explicitly rather than fabricated.
#
# Outputs: tables/anchor_pca_global.csv (sample, PC1-4, group [Set1 subpop label or "3K-RGP reference"])
#          tables/anchor_pca_variance.csv (PC, pct_variance)
COMPENDIUM_ROOT <- Sys.getenv("COMPENDIUM_ROOT", unset = normalizePath(getwd(), winslash = "/"))
source(file.path(COMPENDIUM_ROOT, "analysis/scripts/_paths.R"))
suppressMessages({library(SNPRelate); library(gdsfmt)})

REF_STEM   <- RGP_MERGED_1M
THIN_STEM  <- file.path(INTERIM, "lea", "anchor_global_thin")
GDS_PATH   <- file.path(INTERIM, "anchor_global.gds")

dir.create(dirname(THIN_STEM), showWarnings = FALSE, recursive = TRUE)
cat("=== 3K-RGP global anchor: thin to ~10%, then PCA on the merged 3024-accession panel ===\n")
status <- system2(as.character(PLINK),
                   args = c("--bfile", shQuote(REF_STEM), "--thin", "0.1", "--make-bed",
                             "--seed", "42", "--chr-set", "12", "no-xy", "--allow-extra-chr",
                             "--silent", "--out", shQuote(THIN_STEM)))
if (status != 0 || !file.exists(paste0(THIN_STEM, ".bed")))
  stop(sprintf("PLINK --thin failed (status=%d)", status))
cat(sprintf("thinned to %d SNPs\n", length(readLines(paste0(THIN_STEM, ".bim")))))

if (file.exists(GDS_PATH)) unlink(GDS_PATH)
snpgdsBED2GDS(paste0(THIN_STEM, ".bed"), paste0(THIN_STEM, ".fam"), paste0(THIN_STEM, ".bim"),
              GDS_PATH, verbose = FALSE)
g <- snpgdsOpen(GDS_PATH)
pca <- snpgdsPCA(g, num.thread = 2, verbose = FALSE)
snpgdsClose(g)

vp <- round(pca$varprop[1:6] * 100, 2)
cat(sprintf("n=%d accessions, PC1=%.2f%%  PC2=%.2f%%  PC3=%.2f%%  PC4=%.2f%%\n",
            length(pca$sample.id), vp[1], vp[2], vp[3], vp[4]))

df <- data.frame(sample = pca$sample.id, pca$eigenvect[, 1:4])
names(df)[2:5] <- paste0("PC", 1:4)

# ---- attach Set1 subpopulation labels where available; everything else = reference background ----
subpop <- read.csv(file.path(TAB, "subpop_assignment_set1.csv"), stringsAsFactors = FALSE)
df$group <- "3K-RGP reference (unlabelled)"
m <- match(df$sample, subpop$sample)
df$group[!is.na(m)] <- subpop$subpopulation[m[!is.na(m)]]

n_set1_matched <- sum(!is.na(m))
cat(sprintf("matched %d / %d Set1 accessions to the global PCA by sample ID\n", n_set1_matched, nrow(subpop)))

write.csv(df, file.path(TAB, "anchor_pca_global.csv"), row.names = FALSE)
write.csv(data.frame(PC = paste0("PC", 1:6), pct_variance = vp), file.path(TAB, "anchor_pca_variance.csv"),
          row.names = FALSE)

# ---- quick per-group centroid summary (which region of global diversity each Set1 subpop occupies) ----
agg <- aggregate(cbind(PC1, PC2) ~ group, data = df, FUN = mean)
agg <- agg[order(agg$group), ]
cat("\nPer-group mean PC1/PC2 (reference background vs Set1 subpopulations):\n")
for (i in seq_len(nrow(agg)))
  cat(sprintf("  %-30s PC1=%.4f  PC2=%.4f\n", agg$group[i], agg$PC1[i], agg$PC2[i]))

cat("\n-> tables/anchor_pca_global.csv + anchor_pca_variance.csv\n")
cat("NOTE: Set2 (50K array) NOT anchored -- no pre-merged 3K-RGP+Set2 genotype file found on disk;\n")
cat("      only a common-SNP-position list exists. Documented as a limitation.\n")
