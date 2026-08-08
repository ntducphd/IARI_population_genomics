#!/usr/bin/env python
"""15_fst_wc.py — [stage 15] Weir-Cockerham theta (Fst) among admixture clusters, with
block-bootstrap confidence intervals and a per-cluster sample-size table.

Motivation: Nei's Gst (stage 05) stays as a secondary
estimator; Weir-Cockerham theta is the field-standard, sample-size-corrected estimator and should
lead the differentiation reporting. Implemented with scikit-allel's weir_cockerham_fst (per-locus
variance components a, b, c; global theta = sum(a) / sum(a+b+c), the standard
ratio-of-averages). CIs by block bootstrap over 1-Mb blocks of loci (accounts for LD between
nearby SNPs, which a per-SNP bootstrap would ignore).

Small-cluster rule (reported, not silently applied): clusters with n < 5 are flagged; pairwise
theta involving them is reported but marked low_n = True so the manuscript can state a minimum-n
interpretation rule.

Inputs: {panel}_qc bfile -> VCF via PLINK (--recode vcf-iid), cluster = argmax of the stage-03
admixture Q matrix (same assignment used everywhere).

Outputs (analysis/results/tables/):
  fst_wc_global_{set}.csv    — theta_global, ci95_lo, ci95_hi, n_blocks, n_loci
  fst_wc_pairwise_{set}.csv  — cluster_a, cluster_b, n_a, n_b, theta, ci95_lo, ci95_hi, low_n
  cluster_sizes_{set}.csv    — cluster, n  (feeds the manuscript's minimum-n statement)
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import allel

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PLINK, INTERIM, TAB

PANELS = ["Set1", "Set2"]
BLOCK_BP = 1_000_000
N_BOOT = 200
RNG = np.random.default_rng(42)
LOW_N = 5


def run_plink(*args):
    cmd = [str(PLINK), *[str(a) for a in args], "--chr-set", "12", "no-xy",
           "--allow-extra-chr", "--silent"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("PLINK ERROR:\n", r.stderr[-1500:])
        raise SystemExit(1)


def load_vcf(panel):
    vcf = INTERIM / f"{panel}_qc.vcf"
    if not vcf.exists():
        run_plink("--bfile", INTERIM / f"{panel}_qc", "--recode", "vcf-iid",
                  "--out", INTERIM / f"{panel}_qc")
    cs = allel.read_vcf(str(vcf), fields=["samples", "variants/CHROM", "variants/POS",
                                          "calldata/GT"])
    return (cs["samples"], cs["variants/CHROM"], cs["variants/POS"],
            allel.GenotypeArray(cs["calldata/GT"]))


def cluster_assign(panel, samples):
    q = pd.read_csv(TAB / f"admixture_{panel.lower()}_Q.csv")
    qcols = [c for c in q.columns if c.startswith("Q")]
    m = dict(zip(q["sample"], "C" + q[qcols].to_numpy().argmax(axis=1).astype(str)))
    return np.array([m.get(s, "NA") for s in samples])


def block_ids(chrom, pos):
    return pd.factorize(pd.Series(chrom).astype(str) + "_" +
                        (pd.Series(pos) // BLOCK_BP).astype(str))[0]


def theta_with_ci(gt, subpop_indices, blocks):
    a, b, c = allel.weir_cockerham_fst(gt, subpop_indices)
    num = a.sum(axis=1)                      # per-locus sum of a over alleles
    den = (a + b + c).sum(axis=1)
    keep = np.isfinite(num) & np.isfinite(den)
    num, den, blk = num[keep], den[keep], blocks[keep]
    theta = num.sum() / den.sum()
    # block bootstrap
    uniq = np.unique(blk)
    bnum = np.array([num[blk == u].sum() for u in uniq])
    bden = np.array([den[blk == u].sum() for u in uniq])
    boots = []
    for _ in range(N_BOOT):
        idx = RNG.integers(0, len(uniq), len(uniq))
        boots.append(bnum[idx].sum() / bden[idx].sum())
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(theta), float(lo), float(hi), len(uniq), int(keep.sum())


def main():
    for panel in PANELS:
        samples, chrom, pos, gt = load_vcf(panel)
        clusters = cluster_assign(panel, samples)
        blocks = block_ids(chrom, pos)

        sizes = pd.Series(clusters).value_counts().rename_axis("cluster").reset_index(name="n")
        sizes.to_csv(TAB / f"cluster_sizes_{panel.lower()}.csv", index=False)

        labels = sorted(sizes["cluster"])
        idx = {c: np.where(clusters == c)[0].tolist() for c in labels}

        th, lo, hi, nb, nl = theta_with_ci(gt, [idx[c] for c in labels], blocks)
        pd.DataFrame([{"theta_global": th, "ci95_lo": lo, "ci95_hi": hi,
                       "n_blocks": nb, "n_loci": nl}]).to_csv(
            TAB / f"fst_wc_global_{panel.lower()}.csv", index=False)
        print(f"[15] {panel}: global WC theta = {th:.3f} [{lo:.3f}, {hi:.3f}] "
              f"({nl} loci, {nb} blocks)")

        rows = []
        for i, ca in enumerate(labels):
            for cb in labels[i + 1:]:
                t2, l2, h2, _, _ = theta_with_ci(gt, [idx[ca], idx[cb]], blocks)
                rows.append({"cluster_a": ca, "cluster_b": cb,
                             "n_a": len(idx[ca]), "n_b": len(idx[cb]),
                             "theta": t2, "ci95_lo": l2, "ci95_hi": h2,
                             "low_n": len(idx[ca]) < LOW_N or len(idx[cb]) < LOW_N})
        pd.DataFrame(rows).to_csv(TAB / f"fst_wc_pairwise_{panel.lower()}.csv", index=False)
    print("[15] done")


if __name__ == "__main__":
    main()
