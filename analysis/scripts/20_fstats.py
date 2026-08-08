#!/usr/bin/env python
"""20_fstats.py — [stage 20] Patterson f3 admixture tests among Set 1 subpopulation groups.

Motivation: f-statistics are the standard
human-genetics formalism for admixture (Patterson et al. 2012); applying them to the rice panel
is both a real test (which groups are admixed mixtures of others?) and a cross-disciplinary
bridge the flagship framing calls for. We use the authoritative 3,000 Rice Genomes/SNP-Seek
subpopulation labels (stage 08a) rather than our own clusters where possible — an external,
independent grouping — and fall back to admixture clusters otherwise.

Test: f3(C; A, B) with block-jackknife SE (scikit-allel average_patterson_f3, 1,000-SNP blocks).
Significantly NEGATIVE f3 (Z <= -3) is unambiguous evidence that C is admixed between sources
related to A and B. Positive f3 is not evidence against admixture (drift can mask it) — stated
in the manuscript. Groups with n < 5 are excluded (estimator instability).

Outputs (analysis/results/tables/):
  f3_set1.csv      — target, sourceA, sourceB, f3, se, z, n_blocks, admixed (Z<=-3)
  f3_summary.csv   — groups used with n; count of significant negative tests per target
"""
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import allel

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import INTERIM, TAB

MIN_N = 5
BLEN = 1000


def main():
    cs = allel.read_vcf(str(INTERIM / "Set1_qc.vcf"),
                        fields=["samples", "calldata/GT"])
    samples = cs["samples"]
    gt = allel.GenotypeArray(cs["calldata/GT"])

    sub = pd.read_csv(TAB / "subpop_assignment_set1.csv")
    id_col = sub.columns[0]
    lab_col = [c for c in sub.columns if c.lower() in
               ("subpop", "subpopulation", "label", "group")]
    lab_col = lab_col[0] if lab_col else sub.columns[1]
    m = dict(zip(sub[id_col].astype(str), sub[lab_col].astype(str)))
    labels = np.array([m.get(str(s), "NA") for s in samples])

    counts = pd.Series(labels).value_counts()
    groups = [g for g in counts.index if g != "NA" and counts[g] >= MIN_N]
    print(f"[20] groups used (n>={MIN_N}):",
          ", ".join(f"{g}(n={counts[g]})" for g in groups))

    ac = {g: gt.count_alleles(subpop=np.where(labels == g)[0].tolist()) for g in groups}

    rows = []
    for target in groups:
        for a, b in combinations([g for g in groups if g != target], 2):
            f3, se, z, _, _ = allel.average_patterson_f3(
                ac[target], ac[a], ac[b], blen=BLEN)
            rows.append({"target": target, "sourceA": a, "sourceB": b,
                         "f3": float(f3), "se": float(se), "z": float(z),
                         "admixed": bool(z <= -3)})
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "f3_set1.csv", index=False)

    summ = df.groupby("target")["admixed"].sum().rename("n_sig_negative").reset_index()
    summ["n_tests"] = df.groupby("target").size().values
    summ.to_csv(TAB / "f3_summary.csv", index=False)
    for _, r in summ.iterrows():
        print(f"[20] {r['target']}: {r['n_sig_negative']}/{r['n_tests']} "
              f"tests significantly negative (admixed)")
    print("[20] done -> f3_set1.csv")


if __name__ == "__main__":
    main()
