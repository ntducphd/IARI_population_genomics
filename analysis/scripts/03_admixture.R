# 03_admixture.R — model-based ancestry (sNMF / LEA), K = 1..10, cross-entropy K-selection.
# Deliberately loads ONLY the LEA package (environment audit:
# LEA + adegenet crashed R when loaded together in the same session on this machine, and even some
# LEA/adegenet function introspections crashed non-deterministically. Every pipeline script in this
# compendium therefore loads the minimum package set it actually needs, tested in isolation, and
# is safe to simply re-run if it crashes -- the flakiness observed was not tied to a specific,
# always-reproducible cause).
#
# Pipeline: Set{1,2}_pruned.{bed,bim,fam} (already produced by 01_prep_genotypes.py)
#   -> PLINK --recode vcf (re-derive a VCF from the PRUNED, unlinked SNP set actually used for
#      structure, not the original Set2.vcf from Stage 0 which is pre-prune)
#   -> LEA::vcf2geno
#   -> LEA::snmf, K = 1..10, 5 repetitions, cross-entropy criterion
#   -> best K = argmin mean cross-entropy across repetitions; best run = min cross-entropy at best K
#   -> ancestry (Q) matrix for the best run
#
# Outputs (analysis/results/tables/):
#   admixture_{set}_cv.csv   — K, repetition, cross.entropy (every run, for the K-selection figure)
#   admixture_{set}_Q.csv    — sample + Q1..Q{bestK} ancestry proportions at the best K/run
#   admixture_bestK.csv      — one row per panel: bestK, min mean cross-entropy
COMPENDIUM_ROOT <- Sys.getenv("COMPENDIUM_ROOT", unset = normalizePath(getwd(), winslash = "/"))
source(file.path(COMPENDIUM_ROOT, "analysis/scripts/_paths.R"))
suppressMessages(library(LEA))

KMAX <- 12L; REPS <- 5L; SEED <- 42L
LEADIR <- file.path(INTERIM, "lea"); dir.create(LEADIR, showWarnings = FALSE, recursive = TRUE)

# ---- K-selection: elbow, not naive argmin ---------------------------------------------------
# Cross-entropy in sNMF/ADMIXTURE-style scans on real diversity panels routinely keeps decreasing
# to the edge of the tested range (fine-scale/hierarchical sub-structure keeps improving the fit
# slightly), so "the K with the lowest cross-entropy" is usually just "the largest K tested" and is
# not a meaningful criterion on its own -- confirmed empirically here (see admixture_bestK.csv
# after the first run: both panels picked K=10=KMAX). Report BOTH numbers, honestly:
#   bestK_min   = the K with the global minimum cross-entropy in the tested range (usually KMAX;
#                 reflects fine-scale structure, not "the" answer)
#   bestK_elbow = the smallest K at which the marginal improvement (drop from K to K+1) falls below
#                 `elbow_frac` of the largest single-step drop observed (the conventional ad hoc
#                 elbow/knee rule used throughout the ADMIXTURE/sNMF literature) -- this is the
#                 number to lead with for the coarse "dominant cluster + diverse tail" narrative.
elbow_frac <- 0.10
pick_elbow_K <- function(mean_ce) {
  ord <- mean_ce[order(mean_ce$K), ]
  drops <- -diff(ord$cross_entropy)                      # positive = improvement from K to K+1
  if (length(drops) == 0 || max(drops) <= 0) return(ord$K[1])
  thresh <- elbow_frac * max(drops)
  first_small <- which(drops < thresh)[1]
  if (is.na(first_small)) ord$K[nrow(ord)] else ord$K[first_small]     # K at which the NEXT step is already small
}

panels <- c(Set1 = file.path(INTERIM, "Set1_pruned"), Set2 = file.path(INTERIM, "Set2_pruned"))
bestK_rows <- list()

for (nm in names(panels)) {
  stem <- panels[[nm]]
  cat(sprintf("\n=== [%s] admixture K-scan ===\n", nm))

  # ---- PLINK bed -> vcf (re-derived from the PRUNED set actually used for structure) ----
  vcf_out <- file.path(LEADIR, paste0(nm, "_pruned"))
  status <- system2(as.character(PLINK),
                     args = c("--bfile", shQuote(stem), "--recode", "vcf",
                               "--chr-set", "12", "no-xy", "--allow-extra-chr",
                               "--silent", "--out", shQuote(vcf_out)))
  if (status != 0 || !file.exists(paste0(vcf_out, ".vcf")))
    stop(sprintf("[%s] PLINK --recode vcf failed (status=%d)", nm, status))

  # ---- vcf -> geno ----
  geno_file <- vcf2geno(paste0(vcf_out, ".vcf"), output.file = paste0(vcf_out, ".geno"), force = TRUE)

  # ---- snmf K-scan ----
  proj_dir <- file.path(LEADIR, paste0(nm, "_snmf"))
  if (dir.exists(paste0(vcf_out, ".snmf"))) unlink(paste0(vcf_out, ".snmf"), recursive = TRUE)
  proj <- snmf(geno_file, K = 1:KMAX, repetitions = REPS, project = "new",
               entropy = TRUE, CPU = 2, seed = SEED, iterations = 200, tolerance = 1e-5)

  # ---- cross-entropy per K x repetition ----
  ce_rows <- list(); idx <- 1L
  for (k in 1:KMAX) for (r in 1:REPS) {
    ce_rows[[idx]] <- data.frame(panel = nm, K = k, repetition = r,
                                  cross_entropy = cross.entropy(proj, K = k)[r, 1])
    idx <- idx + 1L
  }
  ce <- do.call(rbind, ce_rows)
  write.csv(ce, file.path(TAB, paste0("admixture_", tolower(nm), "_cv.csv")), row.names = FALSE)

  # ---- K-selection: report both the elbow K (headline) and the global-min K (transparency) ----
  mean_ce <- aggregate(cross_entropy ~ K, data = ce, FUN = mean)
  bestK_elbow <- pick_elbow_K(mean_ce)
  bestK_min   <- mean_ce$K[which.min(mean_ce$cross_entropy)]
  # the Q matrix reported in the main text is built at the elbow K, which is the biologically
  # interpretable "coarse structure" number; bestK_min is kept in admixture_bestK.csv for
  # transparency (as documented: cross-entropy in real diversity panels
  # routinely keeps improving to the edge of the tested range without a sharp global minimum).
  ce_at_bestK <- ce[ce$K == bestK_elbow, ]
  bestrun <- ce_at_bestK$repetition[which.min(ce_at_bestK$cross_entropy)]
  cat(sprintf("[%s] elbow K = %d (mean CE = %.4f) | global-min K = %d (mean CE = %.4f) | best run at elbow K = %d\n",
              nm, bestK_elbow, mean_ce$cross_entropy[mean_ce$K == bestK_elbow],
              bestK_min, min(mean_ce$cross_entropy), bestrun))

  # ---- Q matrix at elbow K / best run ----
  qmat <- Q(proj, K = bestK_elbow, run = bestrun)
  fam <- read.table(paste0(stem, ".fam"), stringsAsFactors = FALSE)
  qdf <- data.frame(sample = fam$V2, qmat)
  names(qdf)[2:(bestK_elbow + 1)] <- paste0("Q", 1:bestK_elbow)
  write.csv(qdf, file.path(TAB, paste0("admixture_", tolower(nm), "_Q.csv")), row.names = FALSE)

  bestK_rows[[nm]] <- data.frame(panel = nm, n = nrow(qdf),
                                  bestK_elbow = bestK_elbow,
                                  mean_ce_at_elbow = round(mean_ce$cross_entropy[mean_ce$K == bestK_elbow], 4),
                                  bestK_global_min = bestK_min,
                                  mean_ce_at_global_min = round(min(mean_ce$cross_entropy), 4))
}

write.csv(do.call(rbind, bestK_rows), file.path(TAB, "admixture_bestK.csv"), row.names = FALSE)
cat("\n-> tables/admixture_{set1,set2}_cv.csv + admixture_{set1,set2}_Q.csv + admixture_bestK.csv\n")
