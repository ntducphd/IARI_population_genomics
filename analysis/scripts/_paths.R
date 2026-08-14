# _paths.R — path config for the R pipeline (mirror of paths.py). `source()` from analysis/scripts/.
# run_all.sh exports COMPENDIUM_ROOT; fallback = current working dir (must be the compendium root).
COMPENDIUM_ROOT <- Sys.getenv("COMPENDIUM_ROOT", unset = normalizePath(getwd(), winslash = "/"))
if (!dir.exists(file.path(COMPENDIUM_ROOT, "analysis")))
  stop("Run from the compendium root, or export COMPENDIUM_ROOT=<compendium root>")

INPUT    <- file.path(COMPENDIUM_ROOT, "analysis/data/input")
INTERIM  <- file.path(COMPENDIUM_ROOT, "analysis/data/interim")
FIG_MAIN <- file.path(COMPENDIUM_ROOT, "analysis/results/figures/main")
FIG_SUPP <- file.path(COMPENDIUM_ROOT, "analysis/results/figures/supp")
TAB      <- file.path(COMPENDIUM_ROOT, "analysis/results/tables")
STAB     <- file.path(COMPENDIUM_ROOT, "analysis/results/supp_tables")
for (d in c(INPUT, INTERIM, FIG_MAIN, FIG_SUPP, TAB, STAB))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---- external data (NOT distributed in this public compendium) ----
# Raw genotype calls, the PLINK 1.9 binary, and the companion phenomic cohort files come from
# restricted/institutional sources outside this repository -- see the manuscript's Data
# Availability statement. Override via env var, or place your own copy under
# analysis/data/external/ (see analysis/data/external/README.md for the expected layout).
EXTERNAL <- Sys.getenv("EXTERNAL_DATA_ROOT", unset = file.path(COMPENDIUM_ROOT, "analysis/data/external"))
.external <- function(env_var, default_rel) {
  v <- Sys.getenv(env_var, unset = NA)
  if (!is.na(v)) return(v)
  file.path(EXTERNAL, default_rel)
}

require_external <- function(path, what) {
  if (!file.exists(path)) {
    stop(sprintf(paste0(
      "[external data missing] %s\n  expected at: %s\n",
      "  This input is not distributed in this public compendium (restricted/institutional ",
      "source -- see the manuscript's Data Availability statement). If you have your own copy, ",
      "place it at that path, or point the matching environment variable at it ",
      "(see analysis/data/external/README.md)."), what, path))
  }
  path
}

PLINK         <- .external("PLINK_BIN", "plink.exe")
GENO          <- .external("RAW_GENOTYPE_DIR", "raw_genotype")
SET1_BED_1M   <- file.path(GENO, "Subset1_150Geno_1M/Genotypes/150genotypes")
SET2_HMP      <- file.path(GENO, "Subset2_147Geno/Genotypes/147SNPgenoypes.hmp.csv")
RGP_MERGED_1M <- file.path(GENO, "Subset1_150Geno_1M/Genotypes/pruned_v2.1")   # 3024 acc x 1.01M (3K-RGP + Set1 merged ref)
.cohort_dir <- .external("PHENOMIC_COHORT_DIR", "phenomic_cohort")
COHORT_SET1 <- file.path(.cohort_dir, "cohort_set1.csv")
COHORT_SET2 <- file.path(.cohort_dir, "cohort_set2.csv")

# GAPIT-Oceanic palette shared across the portfolio (subpopulation / input colours)
PAL_SUBPOP <- c(indica="#2C7FB8", aus="#7FCDBB", `trop-japonica`="#238B45",
                `temp-japonica`="#00441B", aromatic="#D95F0E", admixed="#999999")
