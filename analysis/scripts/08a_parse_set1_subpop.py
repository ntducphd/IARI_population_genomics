#!/usr/bin/env python
"""08a_parse_set1_subpop.py — parse the authoritative 3K-RGP/SNP-Seek subpopulation assignment
for Set1's 150 accessions, already present on disk (mylists-634211649175051595.txt), rather than
re-deriving subpopulation labels from scratch. This is a real, external ground-truth label set
(not our own clustering), which lets Stage 2's admixture/PCA/UMAP clusters be validated AGAINST
known global rice subpopulations rather than only against each other.

Output: tables/subpop_assignment_set1.csv (sample [IRIS_313-xxxx format], subpopulation)
        tables/subpop_capture_set1.csv (subpopulation, n, pct -- the "narrow elite base" evidence)
"""
import re
import sys
from pathlib import Path
from collections import Counter
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import GENO, TAB

SRC = GENO / "Subset1_150Geno_5.2M" / "Wanget al.2018" / "mylists-634211649175051595.txt"

lines = open(SRC, encoding="utf-8").readlines()
rows = []
for line in lines:
    line = line.rstrip("\n")
    if not line.strip() or line.startswith(("VERSION", "VARIETY", "#")):
        continue
    parts = line.split("\t")
    if len(parts) < 8:
        continue
    iris_raw, subpop, country = parts[4], parts[6], parts[7]
    iris_id = iris_raw.strip().replace("IRIS ", "IRIS_")   # "IRIS 313-11155" -> "IRIS_313-11155"
    rows.append(dict(sample=iris_id, subpopulation=subpop.strip(), country=country.strip()))

df = pd.DataFrame(rows)
print(f"parsed {len(df)} Set1 accessions with subpopulation labels")
df.to_csv(TAB / "subpop_assignment_set1.csv", index=False)

counts = df["subpopulation"].value_counts()
capture = counts.reset_index()
capture.columns = ["subpopulation", "n"]
capture["pct"] = round(100 * capture["n"] / capture["n"].sum(), 1)
capture.to_csv(TAB / "subpop_capture_set1.csv", index=False)

print(capture.to_string(index=False))
narrow = capture[capture["subpopulation"].isin(["aus", "indx", "ind2"])]["pct"].sum()
print(f"\naus + indx + ind2 (indica-related groups) = {narrow:.1f}% of Set1 "
      f"-- narrow-elite-base evidence, from authoritative 3K-RGP-derived labels, not our own clustering")
print("\n-> tables/subpop_assignment_set1.csv + subpop_capture_set1.csv")
