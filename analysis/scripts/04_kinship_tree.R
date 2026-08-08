# 04_kinship_tree.R — IBS distance + GRM kinship (SNPRelate) on the LD-pruned genotypes.
# Loads ONLY SNPRelate + gdsfmt (see the environment-audit record environment audit: this R 4.4.3
# build crashes when 3+ compiled-code packages are loaded in one session, even though every PAIR
# tested loads fine -- SNPRelate+gdsfmt+ape crashed, each pair alone did not). The NJ tree itself is
# therefore built in a SEPARATE script, 04b_nj_tree.R, which loads only `ape` and reads the distance
# matrix this script writes -- no joint SNPRelate+ape session anywhere in the pipeline.
#
# Outputs (analysis/results/tables/ + data/interim/):
#   kinship_{set}.csv        — full n x n GRM (GCTA method), long format (sample_a, sample_b, kinship)
#   ibs_dist_{set}.csv       — full n x n IBS-derived distance (1 - IBS), long format, for the tree
#   interim/ibs_dist_{set}.rds — the same distance matrix as an R dist object, for 04b_nj_tree.R
COMPENDIUM_ROOT <- Sys.getenv("COMPENDIUM_ROOT", unset = normalizePath(getwd(), winslash = "/"))
source(file.path(COMPENDIUM_ROOT, "analysis/scripts/_paths.R"))
suppressMessages({library(SNPRelate); library(gdsfmt)})

panels <- c(Set1 = file.path(INTERIM, "Set1_pruned"), Set2 = file.path(INTERIM, "Set2_pruned"))

for (nm in names(panels)) {
  stem <- panels[[nm]]
  gds_path <- paste0(stem, ".gds")
  cat(sprintf("\n=== [%s] kinship + IBS distance ===\n", nm))
  if (!file.exists(gds_path))
    stop(sprintf("[%s] missing %s -- run 02_pca_structure.R first", nm, gds_path))
  g <- snpgdsOpen(gds_path)

  # ---- GRM kinship (GCTA method) -> long format for the heatmap ----
  grm <- snpgdsGRM(g, method = "GCTA", autosome.only = FALSE, verbose = FALSE)
  ids <- grm$sample.id
  K <- grm$grm
  kin_long <- data.frame(
    sample_a = rep(ids, times = length(ids)),
    sample_b = rep(ids, each = length(ids)),
    kinship  = as.vector(K)
  )
  write.csv(kin_long, file.path(TAB, paste0("kinship_", tolower(nm), ".csv")), row.names = FALSE)
  cat(sprintf("[%s] GRM kinship: %d x %d, mean diagonal = %.3f\n", nm, length(ids), length(ids),
              mean(diag(K))))

  # ---- IBS distance (1 - IBS) -> feeds the NJ tree in 04b_nj_tree.R ----
  ibs <- snpgdsIBS(g, autosome.only = FALSE, verbose = FALSE)
  D <- 1 - ibs$ibs
  rownames(D) <- colnames(D) <- ibs$sample.id
  dist_long <- data.frame(
    sample_a = rep(rownames(D), times = ncol(D)),
    sample_b = rep(colnames(D), each = nrow(D)),
    distance = as.vector(D)
  )
  write.csv(dist_long, file.path(TAB, paste0("ibs_dist_", tolower(nm), ".csv")), row.names = FALSE)
  saveRDS(as.dist(D), file.path(INTERIM, paste0("ibs_dist_", tolower(nm), ".rds")))
  cat(sprintf("[%s] IBS distance: %d x %d, mean = %.4f, range = [%.4f, %.4f]\n",
              nm, nrow(D), ncol(D), mean(D[upper.tri(D)]),
              min(D[upper.tri(D)]), max(D[upper.tri(D)])))

  snpgdsClose(g)
}
cat("\n-> tables/kinship_{set}.csv + ibs_dist_{set}.csv + interim/ibs_dist_{set}.rds\n")
