# 02_pca_structure.R — PCA per panel (SNPRelate) on the LD-pruned genotypes.
# Outputs: tables/pca_{set}.csv (sample + PC1..6), tables/pca_variance.csv, and a first quantitative
# read of the "convergent cluster" (dominant-cluster fraction via k-means on the top PCs).
COMPENDIUM_ROOT <- Sys.getenv("COMPENDIUM_ROOT", unset = normalizePath(getwd(), winslash = "/"))
source(file.path(COMPENDIUM_ROOT, "analysis/scripts/_paths.R"))
suppressMessages({library(SNPRelate); library(gdsfmt)})

panels <- c(Set1 = file.path(INTERIM, "Set1_pruned"), Set2 = file.path(INTERIM, "Set2_pruned"))
varrows <- list()
for (nm in names(panels)) {
  stem <- panels[[nm]]; gds <- file.path(INTERIM, paste0(nm, "_pruned.gds"))
  if (file.exists(gds)) unlink(gds)
  snpgdsBED2GDS(paste0(stem, ".bed"), paste0(stem, ".fam"), paste0(stem, ".bim"), gds, verbose = FALSE)
  g <- snpgdsOpen(gds)
  pca <- snpgdsPCA(g, num.thread = 2, verbose = FALSE)
  vp <- round(pca$varprop[1:6] * 100, 2)
  df <- data.frame(sample = pca$sample.id, pca$eigenvect[, 1:6])
  names(df)[2:7] <- paste0("PC", 1:6)
  write.csv(df, file.path(TAB, paste0("pca_", tolower(nm), ".csv")), row.names = FALSE)
  # first read of convergence: k-means (k=3) on top-4 PCs -> size of the dominant cluster
  set.seed(1); km <- kmeans(df[, 2:5], centers = 3, nstart = 25)
  dom <- round(max(table(km$cluster)) / nrow(df) * 100, 1)
  snpgdsClose(g)
  varrows[[nm]] <- data.frame(panel = nm, n = nrow(df), PC1 = vp[1], PC2 = vp[2], PC3 = vp[3],
                              PC4 = vp[4], dom_cluster_pct = dom)
  cat(sprintf("[%s] n=%d  PC1=%.2f%%  PC2=%.2f%%  PC3=%.2f%%  |  dominant k-means cluster = %.1f%% of accessions\n",
              nm, nrow(df), vp[1], vp[2], vp[3], dom))
}
vv <- do.call(rbind, varrows)
write.csv(vv, file.path(TAB, "pca_variance.csv"), row.names = FALSE)
cat("-> tables/pca_{set1,set2}.csv + pca_variance.csv\n")
