# _paths.R — path config for the R pipeline (mirror of paths.py). `source()` from analysis/scripts/.
# run_all.sh exports COMPENDIUM_ROOT; fallback = current working dir (must be the compendium root).
COMPENDIUM_ROOT <- Sys.getenv("COMPENDIUM_ROOT", unset = normalizePath(getwd(), winslash = "/"))
if (!dir.exists(file.path(COMPENDIUM_ROOT, "analysis")))
  stop("Run from the compendium root, or export COMPENDIUM_ROOT=<compendium root>")
WORKSPACE <- normalizePath(file.path(COMPENDIUM_ROOT, ".."), winslash = "/")

INPUT    <- file.path(COMPENDIUM_ROOT, "analysis/data/input")
INTERIM  <- file.path(COMPENDIUM_ROOT, "analysis/data/interim")
FIG_MAIN <- file.path(COMPENDIUM_ROOT, "analysis/results/figures/main")
FIG_SUPP <- file.path(COMPENDIUM_ROOT, "analysis/results/figures/supp")
TAB      <- file.path(COMPENDIUM_ROOT, "analysis/results/tables")
STAB     <- file.path(COMPENDIUM_ROOT, "analysis/results/supp_tables")
for (d in c(INPUT, INTERIM, FIG_MAIN, FIG_SUPP, TAB, STAB))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

PLINK       <- file.path(WORKSPACE, "manuscript_6_gwas_nue/scripts/Plink1.9/plink.exe")
GENO        <- file.path(WORKSPACE, "data/raw/genotype")
SET1_BED_1M <- file.path(GENO, "Subset1_150Geno_1M/Genotypes/150genotypes")
SET2_HMP    <- file.path(GENO, "Subset2_147Geno/Genotypes/147SNPgenoypes.hmp.csv")
COHORT_SET1 <- file.path(WORKSPACE, "manuscript_7_phenomic_selection/analysis/data/input/cohort_set1.csv")
COHORT_SET2 <- file.path(WORKSPACE, "manuscript_7_phenomic_selection/analysis/data/input/cohort_set2.csv")

# GAPIT-Oceanic palette shared across the portfolio (subpopulation / input colours)
PAL_SUBPOP <- c(indica="#2C7FB8", aus="#7FCDBB", `trop-japonica`="#238B45",
                `temp-japonica`="#00441B", aromatic="#D95F0E", admixed="#999999")
