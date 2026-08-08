# 04b_nj_tree.R — neighbour-joining tree from the IBS distance matrix built by 04_kinship_tree.R.
# Loads ONLY `ape` (see 04_kinship_tree.R header: SNPRelate+gdsfmt+ape crashes together on this
# machine, every pair is fine alone -- so the tree-building step is a separate script/process that
# never loads SNPRelate). Bootstrap support (ape::boot.phylo) would require re-deriving IBS
# distances from resampled SNPs, i.e. calling back into SNPRelate from the same session, which is
# exactly the crash-prone combination; bootstrap support is therefore deferred (STRETCH: IQ-TREE
# with ultrafast bootstrap, as documented, if a reviewer requires it).
#
# Outputs:
#   interim/nj_{set}.nwk           — Newick tree
#   tables/nj_{set}_tips.csv       — tip order (for matching figure colouring to subpopulation/admixture)
COMPENDIUM_ROOT <- Sys.getenv("COMPENDIUM_ROOT", unset = normalizePath(getwd(), winslash = "/"))
source(file.path(COMPENDIUM_ROOT, "analysis/scripts/_paths.R"))
suppressMessages(library(ape))

panels <- c("Set1", "Set2")
for (nm in panels) {
  cat(sprintf("\n=== [%s] NJ tree ===\n", nm))
  dist_path <- file.path(INTERIM, paste0("ibs_dist_", tolower(nm), ".rds"))
  if (!file.exists(dist_path))
    stop(sprintf("[%s] missing %s -- run 04_kinship_tree.R first", nm, dist_path))
  D <- readRDS(dist_path)

  tr <- nj(D)
  tr <- ladderize(tr)

  nwk_path <- file.path(INTERIM, paste0("nj_", tolower(nm), ".nwk"))
  write.tree(tr, file = nwk_path)

  tips_path <- file.path(TAB, paste0("nj_", tolower(nm), "_tips.csv"))
  write.csv(data.frame(tip_order = seq_along(tr$tip.label), sample = tr$tip.label),
            tips_path, row.names = FALSE)

  cat(sprintf("[%s] NJ tree: %d tips, total branch length = %.4f\n",
              nm, length(tr$tip.label), sum(tr$edge.length)))
  cat(sprintf("  -> %s\n  -> %s\n", nwk_path, tips_path))
}
cat("\n-> interim/nj_{set}.nwk + tables/nj_{set}_tips.csv\n")
