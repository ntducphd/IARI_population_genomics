# External data (not shipped in this compendium)

This directory is intentionally empty in the public compendium. Several pipeline stages need
inputs from restricted/institutional sources or a companion compendium — see the manuscript's
Data Availability statement for access details. `analysis/scripts/paths.py` (Python) and
`_paths.R` (R) resolve these paths and raise a clear `[external data missing]` error, naming the
exact expected path, if a stage runs without them.

Either place your own copy at the path shown below, or set the matching environment variable to
point anywhere on your machine.

| Env var | Default path (under this folder) | What it is |
|---|---|---|
| `PLINK_BIN` | `plink.exe` | PLINK 1.9 binary |
| `RAW_GENOTYPE_DIR` | `raw_genotype/` | Raw genotype calls (Set 1 WGS bfile, Set 2 50K HapMap, 3K-RGP reference genotypes) |
| `PHENOMIC_COHORT_DIR` | `phenomic_cohort/` | `cohort_set1.csv` / `cohort_set2.csv` from the companion phenomic-selection compendium (public — see the manuscript for the citation) |

You can also override the whole tree at once with `EXTERNAL_DATA_ROOT`.
