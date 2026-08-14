#!/usr/bin/env bash
# restore_env.sh — rebuild the analysis environment (3 phases).
# Usage: bash restore_env.sh    (requires conda; then `conda activate population-genomics-compendium`)
set -e

echo "[1/3] conda environment (Python + R + Java) from environment.yml"
conda env create -f environment.yml || conda env update -f environment.yml

echo "[2/3] R packages (CRAN, version-pinned via remotes)"
conda run -n population-genomics-compendium Rscript - <<'RS'
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
pins <- c(rehh = "3.2.3", pcadapt = "4.4.1", adegenet = "2.1.11",
          pegas = "1.4", ape = "5.8.1", vcfR = "1.16.0")
for (p in names(pins)) {
  if (!requireNamespace(p, quietly = TRUE) ||
      as.character(packageVersion(p)) != pins[[p]]) {
    remotes::install_version(p, version = pins[[p]], upgrade = "never")
  }
}
RS

echo "[3/3] Bioconductor packages (needed only to re-run stages 01-12 from raw genotypes)"
conda run -n population-genomics-compendium Rscript - <<'RS'
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SNPRelate", "gdsfmt", "LEA"), update = FALSE, ask = FALSE)
RS

echo "Done. Activate with: conda activate population-genomics-compendium"
echo "Third-party jars (Beagle 5.5, Stairway Plot v2.1.2) are fetched to analysis/tools/ by"
echo "their stage scripts; PLINK 1.9 is an external, restricted-access install -- see"
echo "analysis/data/external/README.md (PLINK_BIN env var)."
