# 05b_amova.R — hierarchical AMOVA (pegas) on the IBS distance from Stage 3, grouped by the
# admixture-based cluster from Stage 2. Loads ONLY `pegas` (see 04_kinship_tree.R / 05_fst_amova.py
# headers for the established rule on this machine: keep each script's package set minimal and
# tested). hierfstat/poppr were tried for this step and abandoned -- hierfstat crashed R
# reproducibly on trivial operations (documented environment audit); poppr's amova path goes
# through adegenet genind objects, the same family that crashed on read.PLINK. pegas::amova()
# works directly on a distance matrix + a grouping data.frame, with no adegenet dependency.
#
# Group labels: same admixture argmax-Q elbow-K cluster used in 05_fst_amova.py's Fst calculation
# (data-driven proxy for subpopulation until Stage 7's formal 3K-RGP XI/GJ/cA/cB assignment).
#
# Outputs: tables/amova_{set}.csv — variance components (among-group, within-group) as % of total,
#          plus the permutation-based significance test pegas::amova() runs internally.
COMPENDIUM_ROOT <- Sys.getenv("COMPENDIUM_ROOT", unset = normalizePath(getwd(), winslash = "/"))
source(file.path(COMPENDIUM_ROOT, "analysis/scripts/_paths.R"))
suppressMessages(library(pegas))

panels <- c("Set1", "Set2")
NPERM <- 999L

rows <- list()
for (nm in panels) {
  cat(sprintf("\n=== [%s] AMOVA ===\n", nm))
  dist_path <- file.path(INTERIM, paste0("ibs_dist_", tolower(nm), ".rds"))
  q_path <- file.path(TAB, paste0("admixture_", tolower(nm), "_Q.csv"))
  if (!file.exists(dist_path)) stop(sprintf("[%s] missing %s -- run 04_kinship_tree.R first", nm, dist_path))
  if (!file.exists(q_path)) stop(sprintf("[%s] missing %s -- run 03_admixture.R first", nm, q_path))

  D <- readRDS(dist_path)
  qdf <- read.csv(q_path, stringsAsFactors = FALSE)
  qcols <- grep("^Q", names(qdf), value = TRUE)
  cluster <- apply(as.matrix(qdf[, qcols]), 1, which.max)
  names(cluster) <- qdf$sample

  # align the distance matrix's label order to the cluster vector
  dlabs <- labels(D)
  cluster <- cluster[dlabs]
  if (any(is.na(cluster))) stop(sprintf("[%s] cluster labels do not fully match the distance matrix's samples", nm))

  strata_df <- data.frame(cluster = factor(cluster))
  set.seed(42)
  res <- amova(D ~ cluster, data = strata_df, nperm = NPERM, is.squared = FALSE)

  # pegas amova returns $tab with SSD, MSD, df per stratum (Between-samples / Within-samples), and
  # $varcomp with the variance components + their %. Extract both robustly.
  varcomp <- res$varcomp
  pct <- round(100 * varcomp$sigma2 / sum(varcomp$sigma2), 2)
  cat(sprintf("[%s] variance components: %s\n", nm,
              paste(sprintf("%s=%.2f%%", rownames(varcomp), pct), collapse = ", ")))

  rows[[nm]] <- data.frame(panel = nm, n = length(dlabs), n_clusters = nlevels(strata_df$cluster),
                            stratum = rownames(varcomp), sigma2 = round(varcomp$sigma2, 5),
                            pct_variance = pct)
  write.csv(rows[[nm]], file.path(TAB, paste0("amova_", tolower(nm), ".csv")), row.names = FALSE)
}
cat("\n-> tables/amova_{set}.csv\n")
