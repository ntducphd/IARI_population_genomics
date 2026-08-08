#!/usr/bin/env python
"""12_audit.py — QC gate for this compendium. Checks that every figure/table the MANIFEST
promises actually exists (PNG+PDF pairs for figures), and runs basic sanity checks on the key
numbers (p-values in [0,1], correlations in [-1,1], no unexpected NaNs in headline results).
Exits 1 if any FAIL is found (mirrors the house convention, e.g. manuscript_6's check_consistency.py).
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import TAB, STAB, FIG_MAIN, FIG_SUPP, INTERIM

PASSED, WARNINGS, ERRORS = [], [], []
def ok(m): PASSED.append(m)
def warn(m): WARNINGS.append(m)
def err(m): ERRORS.append(m)

# ---- 1. figure existence (PNG + PDF pairs) ----
EXPECTED_MAIN = ["Fig01_design", "Fig02_structure", "Fig03_tree_kinship",
                 "Fig04_diff_diversity_ld", "Fig05_concordance_mechanism",
                 "Fig06_genome_history", "Fig07_selection", "Fig08_application",
                 "Fig09_mechanism_visible", "Fig10_platforms"]
EXPECTED_SUPP = ["SuppFig01_umap_consensus", "SuppFig02_anchor_3krgp", "SuppFig03_stairway",
                 "SuppFig04_f3", "SuppFig05_dapc"]
for stem in EXPECTED_MAIN:
    png, pdf = FIG_MAIN / f"{stem}.png", FIG_MAIN / f"{stem}.pdf"
    if png.exists() and pdf.exists():
        ok(f"Fig {stem}: PNG+PDF present")
    else:
        err(f"Fig {stem}: MISSING ({'no PNG' if not png.exists() else ''} {'no PDF' if not pdf.exists() else ''})")
for stem in EXPECTED_SUPP:
    png, pdf = FIG_SUPP / f"{stem}.png", FIG_SUPP / f"{stem}.pdf"
    if png.exists() and pdf.exists():
        ok(f"{stem}: PNG+PDF present")
    else:
        err(f"{stem}: MISSING")

# ---- 2. table existence ----
EXPECTED_TABLES = ["Table_1_panel_diversity_summary", "Table_2_fst_amova", "Table_3_ld_decay",
                   "Table_4_pillarB_concordance", "Table_5_core_collection",
                   "Table_6_genome_history", "Table_7_selection"]
for stem in EXPECTED_TABLES:
    md = TAB / f"{stem}.md"
    if md.exists():
        ok(f"Table {stem}.md present")
    else:
        err(f"Table {stem}.md MISSING")

EXPECTED_SUPP_TABLES = [f"SuppTable_S{i}_{n}" for i, n in enumerate(
    ["prep_summary", "subpop_capture_set1", "pairwise_fst", "maf_spectrum",
     "structure_consensus", "core_collection_curve", "underused_accessions",
     "f3_admixture_tests", "allelic_richness", "roh_by_cluster",
     "robustness_summary", "fdr_ledger", "confound_tests", "classifier_robustness",
     "ne_interval", "overlap_null_validation", "pst_fst_sweep", "temporal_mantel",
     "classifier_uncertainty", "mechanism_partition"], start=1)]
for stem in EXPECTED_SUPP_TABLES:
    csv = STAB / f"{stem}.csv"
    if csv.exists():
        ok(f"{stem}.csv present")
    else:
        err(f"{stem}.csv MISSING")

# ---- 3. sanity checks on headline numbers ----
def check_p_range(path, col, label):
    if not path.exists():
        warn(f"{label}: file missing, skip p-value range check"); return
    df = pd.read_csv(path)
    if col not in df.columns:
        warn(f"{label}: column {col} missing"); return
    bad = df[(df[col] < 0) | (df[col] > 1)]
    if len(bad):
        err(f"{label}: {len(bad)} rows with {col} outside [0,1]")
    else:
        ok(f"{label}: all {col} values in [0,1] ({len(df)} rows)")

for panel in ["set1", "set2"]:
    check_p_range(TAB / f"concordance_mantel_{panel}.csv", "p", f"Mantel p ({panel})")
    check_p_range(TAB / f"fst_global_{panel}.csv", "perm_p", f"Fst permutation p ({panel})")
    check_p_range(TAB / f"concordance_procrustes_{panel}.csv", "perm_p", f"Procrustes p ({panel})")

def check_corr_range(path, col, label):
    if not path.exists():
        warn(f"{label}: file missing"); return
    df = pd.read_csv(path)
    bad = df[(df[col] < -1) | (df[col] > 1)]
    if len(bad):
        err(f"{label}: {len(bad)} rows with {col} outside [-1,1]")
    else:
        ok(f"{label}: all {col} values in [-1,1] ({len(df)} rows)")

for panel in ["set1", "set2"]:
    check_corr_range(TAB / f"concordance_mantel_{panel}.csv", "r", f"Mantel r ({panel})")

# ---- 4. headline result cross-check: is genomic<->phenomic Mantel significant in BOTH panels? ----
mant1 = pd.read_csv(TAB / "concordance_mantel_set1.csv")
mant2 = pd.read_csv(TAB / "concordance_mantel_set2.csv")
r1 = mant1[mant1.comparison == "genomic~phenomic"].iloc[0]
r2 = mant2[mant2.comparison == "genomic~phenomic"].iloc[0]
if r1.p < 0.05 and r2.p < 0.05:
    ok(f"PILLAR B headline result replicates across both disjoint panels: "
       f"Set1 r={r1.r:.3f} p={r1.p:.3f}, Set2 r={r2.r:.3f} p={r2.p:.3f}")
else:
    err(f"PILLAR B headline result does NOT replicate in both panels (Set1 p={r1.p}, Set2 p={r2.p})")

# ---- 5. no NaN in key summary tables ----
for stem in ["Table_1_panel_diversity_summary", "Table_2_fst_amova", "Table_5_core_collection"]:
    df = pd.read_csv(TAB / f"{stem}.csv")
    n_nan = df.isna().sum().sum()
    if n_nan == 0:
        ok(f"{stem}.csv: no unexpected NaN")
    else:
        warn(f"{stem}.csv: {n_nan} NaN cell(s) -- verify expected (e.g. Ne_LD_based note)")

# ---- 6. Set2 3K-RGP anchor gap explicitly flagged (not silently missing) ----
if not (TAB / "subpop_assignment_set2.csv").exists():
    warn("Set2 3K-RGP subpopulation assignment: confirmed absent, as documented "
         "(a deliberate, documented gap)")
else:
    warn("Set2 3K-RGP subpopulation assignment file now exists -- update "
         "this check if the gap has been closed")

# ---- report ----
print("=" * 70)
print("COMPENDIUM QC AUDIT")
print("=" * 70)
print(f"\nPASSED ({len(PASSED)}):")
for m in PASSED:
    print(f"  OK  {m}")
if WARNINGS:
    print(f"\nWARNINGS ({len(WARNINGS)}):")
    for m in WARNINGS:
        print(f"  WARN  {m}")
if ERRORS:
    print(f"\nFAILED ({len(ERRORS)}):")
    for m in ERRORS:
        print(f"  FAIL  {m}")
    print(f"\nAudit FAILED: {len(ERRORS)} error(s), {len(WARNINGS)} warning(s), {len(PASSED)} passed.")
    sys.exit(1)
else:
    print(f"\nAudit PASSED: {len(PASSED)} checks, {len(WARNINGS)} warnings, 0 errors.")
    sys.exit(0)
