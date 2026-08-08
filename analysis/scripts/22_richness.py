#!/usr/bin/env python
"""22_richness.py — [stage 22] rarefied allelic richness and private alleles per admixture
cluster, both panels.

Motivation: allelic richness and private alleles are
the expected companions of He/Ho in any germplasm diversity table; richness must be
rarefaction-corrected because clusters differ in size (larger samples find more alleles
mechanically).

Method: per biallelic locus and cluster, expected number of distinct alleles in a rarefied sample
of g gene copies (Hurlbert 1971; Kalinowski 2004):
    A_g = sum_over_alleles [ 1 - C(N - N_i, g) / C(N, g) ]
with g = 2 * min cluster n (gene copies), averaged over loci -> mean rarefied richness in [1, 2].
Private alleles: count of alleles observed in exactly one cluster (on the QC'd SNP set),
attributed to that cluster; reported raw and per-1,000 SNPs.

Outputs (analysis/results/tables/):
  richness_{set}.csv  — cluster, n, mean_rarefied_richness_Ag, n_private_alleles,
                         private_per_1k_snps
  richness_summary.csv — per panel: g used, n_clusters, min/max richness, top private cluster
"""
import sys
from math import lgamma
from pathlib import Path

import numpy as np
import pandas as pd
import allel

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import INTERIM, TAB

MIN_N = 3


def log_choose(n, k):
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def rarefied_richness(ac, g):
    """Mean expected #alleles at rarefaction depth g gene copies, over loci with N >= g."""
    N = ac.sum(axis=1)
    keep = N >= g
    acn, Nn = ac[keep], N[keep]
    total = np.zeros(len(Nn), dtype=float)
    for a in range(acn.shape[1]):
        Ni = acn[:, a]
        with np.errstate(all="ignore"):
            lc = np.array([log_choose(nv - ni, g) - log_choose(nv, g)
                           if nv - ni >= g else -np.inf
                           for nv, ni in zip(Nn, Ni)])
        total += 1.0 - np.exp(lc)
    return float(total.mean()), int(keep.sum())


def main():
    summary = []
    for panel in ["Set1", "Set2"]:
        cs = allel.read_vcf(str(INTERIM / f"{panel}_qc.vcf"),
                            fields=["samples", "calldata/GT"])
        samples = cs["samples"]
        gt = allel.GenotypeArray(cs["calldata/GT"])

        q = pd.read_csv(TAB / f"admixture_{panel.lower()}_Q.csv")
        qcols = [c for c in q.columns if c.startswith("Q")]
        m = dict(zip(q["sample"], "C" + q[qcols].to_numpy().argmax(axis=1).astype(str)))
        labels = np.array([m.get(s, "NA") for s in samples])
        counts = pd.Series(labels).value_counts()
        clusters = [c for c in sorted(counts.index) if c != "NA" and counts[c] >= MIN_N]

        g = 2 * int(min(counts[c] for c in clusters))
        ac = {c: gt.count_alleles(subpop=np.where(labels == c)[0].tolist())
              for c in clusters}

        # private alleles: allele present (count>0) in exactly one cluster
        present = np.stack([ac[c] > 0 for c in clusters])           # (k, loci, alleles)
        n_present = present.sum(axis=0)
        rows = []
        for i, c in enumerate(clusters):
            priv = int((present[i] & (n_present == 1)).sum())
            ag, n_loci = rarefied_richness(np.asarray(ac[c]), g)
            rows.append({"cluster": c, "n": int(counts[c]),
                         "mean_rarefied_richness_Ag": ag,
                         "n_private_alleles": priv,
                         "private_per_1k_snps": 1000.0 * priv / gt.shape[0]})
        df = pd.DataFrame(rows)
        df.to_csv(TAB / f"richness_{panel.lower()}.csv", index=False)
        top = df.loc[df["n_private_alleles"].idxmax()]
        summary.append({"panel": panel, "g_gene_copies": g, "n_clusters": len(clusters),
                        "richness_min": df["mean_rarefied_richness_Ag"].min(),
                        "richness_max": df["mean_rarefied_richness_Ag"].max(),
                        "top_private_cluster": top["cluster"],
                        "top_private_n": int(top["n_private_alleles"])})
        print(f"[22] {panel}: g={g}, richness {df['mean_rarefied_richness_Ag'].min():.3f}-"
              f"{df['mean_rarefied_richness_Ag'].max():.3f}, "
              f"most private alleles: {top['cluster']} ({int(top['n_private_alleles'])})")
    pd.DataFrame(summary).to_csv(TAB / "richness_summary.csv", index=False)
    print("[22] done -> richness_summary.csv")


if __name__ == "__main__":
    main()
