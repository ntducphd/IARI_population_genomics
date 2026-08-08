#!/usr/bin/env python
# paths.py — single source of truth for every path in this compendium.
# Self-locating: works no matter where the compendium is moved (rename/move-proof).
# Import in any pipeline script:  from paths import *
from pathlib import Path
import sys

ROOT      = Path(__file__).resolve().parents[2]          # manuscript_12_population_genomics/
WORKSPACE = ROOT.parent                                   # the parent workspace folder
SCRIPTS   = ROOT / "analysis/scripts"

# ---- compendium data/results ----
INPUT     = ROOT / "analysis/data/input"                  # harmonised, committed inputs
INTERIM   = ROOT / "analysis/data/interim"                # PLINK/GDS intermediates (regenerable)
FIG       = ROOT / "analysis/results/figures"
FIG_MAIN  = FIG / "main";  FIG_SUPP = FIG / "supp"
TAB       = ROOT / "analysis/results/tables"
STAB      = ROOT / "analysis/results/supp_tables"
for _d in (INPUT, INTERIM, FIG_MAIN, FIG_SUPP, TAB, STAB): _d.mkdir(parents=True, exist_ok=True)

# ---- tools (verified 2026-08-03; `where` returns not-found -> use full paths) ----
# 2026-08-08: the compendium runs on more than one machine (different R versions / user profiles),
# so RSCRIPT falls back from the documented 4.4.3 path to whatever Rscript is installed/on PATH.
PLINK   = WORKSPACE / "manuscript_6_gwas_nue/scripts/Plink1.9/plink.exe"


def _find_rscript():
    import shutil
    cands = [Path("C:/Program Files/R/R-4.4.3/bin/Rscript.exe")]
    pf = Path("C:/Program Files/R")
    if pf.exists():
        cands += sorted(pf.glob("R-*/bin/Rscript.exe"), reverse=True)
    for c in cands:
        if c.exists():
            return c
    w = shutil.which("Rscript")
    return Path(w) if w else cands[0]


RSCRIPT = _find_rscript()
PYTHON  = Path(sys.executable) if sys.executable else WORKSPACE / ".venv/Scripts/python.exe"

# ---- raw genotype sources (read-only, in the shared data lake) ----
GENO = WORKSPACE / "data/raw/genotype"
SET1_BED_1M   = GENO / "Subset1_150Geno_1M/Genotypes/150genotypes"          # 150 acc x 1.01M SNP (PLINK bfile stem)
SET1_BED_5M   = GENO / "Subset1_150Geno_5.2M/Genotypes/VTrice"              # 150 acc x 5.23M (dense; fine-scale)
SET2_HMP      = GENO / "Subset2_147Geno/Genotypes/147SNPgenoypes.hmp.csv"   # 147 acc x ~50K (HapMap)
RGP_SNPLIST   = GENO / "Subset1_150Geno_5.2M/Wanget al.2018/3K-HDRA-snp-comm-miss5pc.txt"  # HDRA common-SNP positions
RGP_MERGED_1M = GENO / "Subset1_150Geno_1M/Genotypes/pruned_v2.1"           # 3024 acc x 1.01M (3K-RGP + Set1 merged ref)

# ---- cross-paper phenotype/phenomic (Pillar B) — cohort files already ID-aligned per panel ----
P7 = WORKSPACE / "manuscript_7_phenomic_selection/analysis/data/input"
COHORT_SET1 = P7 / "cohort_set1.csv"     # Set1: genotype IDs + traditional traits + 204 phenomic features
COHORT_SET2 = P7 / "cohort_set2.csv"     # Set2: idem

PANELS = {
    "Set1": {"geno": SET1_BED_1M, "kind": "bfile", "cohort": COHORT_SET1, "platform": "WGS (1.01M)"},
    "Set2": {"geno": SET2_HMP,    "kind": "hmp",   "cohort": COHORT_SET2, "platform": "50K array"},
}

if __name__ == "__main__":
    print("ROOT      :", ROOT)
    for k, v in {"PLINK": PLINK, "RSCRIPT": RSCRIPT, "PYTHON": PYTHON,
                 "SET1_BED_1M.bed": SET1_BED_1M.with_suffix(".bed"),
                 "SET2_HMP": SET2_HMP, "RGP_SNPLIST": RGP_SNPLIST,
                 "COHORT_SET1": COHORT_SET1, "COHORT_SET2": COHORT_SET2}.items():
        print(f"  {'OK ' if Path(str(v)).exists() else 'MISS'} {k}: {v}")
