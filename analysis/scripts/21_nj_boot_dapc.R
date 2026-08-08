#!/usr/bin/env Rscript
# 21_nj_boot_dapc.R -- [stage 21] (a) bootstrap support for the neighbour-joining trees;
# (b) DAPC retry as the fourth structure-inference method.
#
# Motivation: (NJ trees shipped without support values -> no interior-branch
# claim is currently defensible) and T2-8 (the earlier adegenet/DAPC "deterministic crash" has the
# signature of the documented Git-Bash/MSYS2 SVD segfault on this Windows setup, not an adegenet
# bug -- this script re-attempts it under PowerShell; on success DAPC joins the cross-method
# consensus, on failure the documented drop stands with the retry recorded).
#
# (a) NJ bootstrap: 100 SNP-resampling replicates on the LD-pruned dosage matrix; per replicate,
#     Euclidean-dosage distance -> ape::nj; support = ape::prop.clades on the original tree.
#     (Euclidean on dosages is monotone-equivalent to 1-IBS for biallelic SNPs up to missingness;
#     the original tree is rebuilt on the same matrix so support values refer to its own topology.)
# (b) DAPC: adegenet::df2genind is bypassed -- we feed dosages directly to dapc.data.frame via
#     find.clusters (K fixed at each panel's elbow K from stage 03), n.pca chosen by
#     a-score-free heuristic n/3 (adegenet manual guidance), n.da = K-1.
#
# Outputs (analysis/results/tables/):
#   nj_boot_support_{set}.csv — node, bootstrap_support (0-100)
#   nj_boot_summary.csv       — per panel: n_nodes, median/share>=70 support
#   dapc_assign_{set}.csv     — sample, dapc_cluster   (only if DAPC succeeds)
#   dapc_status.csv           — panel, status (ok / failed: <message>)
# Run via PowerShell only (SVD house rule).

suppressMessages({ library(ape) })

args      <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPTS   <- dirname(normalizePath(this_file))
ROOT      <- normalizePath(file.path(SCRIPTS, "..", ".."))
INTERIM   <- file.path(ROOT, "analysis", "data", "interim")
TAB       <- file.path(ROOT, "analysis", "results", "tables")
PLINK     <- file.path(dirname(ROOT), "manuscript_6_gwas_nue", "scripts", "Plink1.9", "plink.exe")

ELBOW_K <- c(Set1 = 7L, Set2 = 9L)
N_BOOT  <- 100L
set.seed(42)

read_dosage <- function(panel) {
  raw <- file.path(INTERIM, paste0(panel, "_pruned.raw"))
  if (!file.exists(raw)) {
    status <- system2(PLINK, c("--bfile", file.path(INTERIM, paste0(panel, "_pruned")),
                               "--recode", "A", "--chr-set", "12", "no-xy",
                               "--allow-extra-chr", "--silent",
                               "--out", file.path(INTERIM, paste0(panel, "_pruned"))))
    if (status != 0) stop("PLINK recode A failed for ", panel)
  }
  d <- read.table(raw, header = TRUE, check.names = FALSE)
  ids <- d$IID
  x <- as.matrix(d[, -(1:6)])
  # mean-impute missing dosages per SNP (same policy as the stage-03b UMAP input)
  cm <- colMeans(x, na.rm = TRUE)
  idx <- which(is.na(x), arr.ind = TRUE)
  if (nrow(idx)) x[idx] <- cm[idx[, 2]]
  rownames(x) <- ids
  x
}

boot_rows <- list(); dapc_status <- list()

for (panel in c("Set1", "Set2")) {
  x <- read_dosage(panel)

  ## (a) NJ bootstrap ---------------------------------------------------------
  tr0 <- nj(dist(x))
  boots <- vector("list", N_BOOT)
  for (b in seq_len(N_BOOT)) {
    cols <- sample.int(ncol(x), replace = TRUE)
    boots[[b]] <- nj(dist(x[, cols]))
  }
  supp <- prop.clades(tr0, boots, rooted = FALSE)
  supp[is.na(supp)] <- 0
  supp_pct <- round(100 * supp / N_BOOT)
  write.csv(data.frame(node = seq_along(supp_pct) + Ntip(tr0),
                       bootstrap_support = supp_pct),
            file.path(TAB, paste0("nj_boot_support_", tolower(panel), ".csv")),
            row.names = FALSE)
  write.tree(tr0, file.path(INTERIM, paste0("nj_", tolower(panel), "_boot.nwk")))
  boot_rows[[panel]] <- data.frame(panel = panel, n_internal_nodes = length(supp_pct),
                                   median_support = median(supp_pct),
                                   share_ge70 = mean(supp_pct >= 70))
  cat(sprintf("[21] %s NJ bootstrap: median support %d%%, %.0f%% of nodes >= 70%%\n",
              panel, median(supp_pct), 100 * mean(supp_pct >= 70)))

  ## (b) DAPC retry -----------------------------------------------------------
  st <- tryCatch({
    suppressMessages(library(adegenet))
    K <- ELBOW_K[[panel]]
    grp <- find.clusters(x, n.pca = min(50, nrow(x) - 1), n.clust = K)
    fit <- dapc(x, grp$grp, n.pca = round(nrow(x) / 3), n.da = K - 1)
    write.csv(data.frame(sample = rownames(x),
                         dapc_cluster = paste0("D", as.integer(fit$grp) - 1L)),
              file.path(TAB, paste0("dapc_assign_", tolower(panel), ".csv")),
              row.names = FALSE)
    cat(sprintf("[21] %s DAPC: SUCCESS (K = %d) -- the earlier crash was environmental\n",
                panel, K))
    "ok"
  }, error = function(e) paste("failed:", conditionMessage(e)))
  dapc_status[[panel]] <- data.frame(panel = panel, status = st)
}

write.csv(do.call(rbind, boot_rows), file.path(TAB, "nj_boot_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, dapc_status), file.path(TAB, "dapc_status.csv"), row.names = FALSE)
cat("[21] done\n")
