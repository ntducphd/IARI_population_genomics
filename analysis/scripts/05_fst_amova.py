#!/usr/bin/env python
"""05_fst_amova.py — pairwise and global Fst between admixture-based clusters (Nei's Gst,
Hs/Ht formulation), computed directly in Python.

Why not hierfstat (R)? hierfstat is installed and loads fine alone, but repeatedly crashed R
(reproducibly, across multiple retries) the moment its bundled example dataset was inspected
(even a plain write.csv() after data(gtrunchier) crashed).
Reimplementing a textbook Fst formula directly, in a language/environment with no such fragility,
is more reliable than fighting an unstable R package's undocumented-here input format. Nei's Gst
(1973/1987) is the classic, easily-verified Hs/Ht estimator:
    Ht (locus) = 2 * pbar * (1 - pbar),   pbar = pooled/weighted allele frequency across groups
    Hs (locus) = mean over groups of 2 * p_g * (1 - p_g)
    Fst (locus) = 1 - Hs / Ht
Genome-wide/pairwise Fst uses the ratio-of-averages form (mean Hs / mean Ht across loci), which is
more stable than averaging per-locus Fst directly (avoids division blow-ups at near-monomorphic
loci) and is the standard way multi-locus Fst is reported.

Group labels: the admixture argmax-Q cluster from 03_admixture.R / 03c_structure_consensus.py (the
elbow-K clustering already validated by ARI/NMI against PCA and UMAP in Stage 2). This is a
data-driven proxy for subpopulation until the formal 3K-RGP XI/GJ/cA/cB assignment (Stage 7) is
built; Stage 7's labels should be cross-referenced against this grouping once available, not
treated as a second, independent labeling to reconcile from scratch.

Permutation test: group labels are shuffled (keeping the genotype matrix fixed) `n_perm` times to
build a null distribution for the global (multi-group) Fst, giving a permutation p-value without
needing hierfstat/StAMPP's bootstrap machinery.

Outputs (analysis/results/tables/):
  fst_pairwise_{set}.csv   — cluster_a, cluster_b, n_a, n_b, fst
  fst_global_{set}.csv     — n_clusters, global_fst, perm_p (n_perm permutations)
"""
import subprocess
import sys
import warnings
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PLINK, INTERIM, TAB

# Benign, expected: a small admixture cluster can be 100% missing at a handful of loci, giving an
# all-NaN slice for that (locus, group) combination. nanmean warns per slice; the result is NaN,
# which the `ok` mask below correctly excludes from the Hs/Ht ratio -- not a silent error.
warnings.filterwarnings("ignore", message="Mean of empty slice")

SEED = 42
N_PERM = 999
PANELS = ["Set1", "Set2"]
rng = np.random.default_rng(SEED)


def load_dosage(panel):
    # LD-pruned set (matches PCA/admixture/kinship in Stage 2-3): using the full QC set here would
    # both be ~20x slower for the permutation test AND statistically questionable, since physically
    # linked SNPs are not independent draws and would understate the true sampling variance in the
    # permutation null (the same reason PCA/admixture/kinship all use the pruned set, not QC).
    stem = INTERIM / f"{panel}_pruned"
    out = INTERIM / "lea" / f"{panel}_fst_raw"
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [str(PLINK), "--bfile", str(stem), "--recode", "A",
           "--chr-set", "12", "no-xy", "--allow-extra-chr", "--silent", "--out", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("PLINK ERROR:\n", r.stderr[-1500:])
        raise SystemExit(1)
    df = pd.read_csv(out.with_suffix(".raw"), sep=r"\s+")
    meta = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]
    snp_cols = [c for c in df.columns if c not in meta]
    X = df[snp_cols].to_numpy(dtype=float)   # 0/1/2 dosage, NaN = missing
    return df["IID"].to_numpy(), X


def group_hs_ht(X, labels):
    """Per-locus Ht (pooled) and Hs (mean within-group), ignoring missing per column."""
    groups = np.unique(labels)
    n_loci = X.shape[1]
    # pooled allele frequency per locus (mean dosage / 2), ignoring NaN
    pbar = np.nanmean(X, axis=0) / 2.0
    ht = 2 * pbar * (1 - pbar)

    hs_per_group = np.zeros((len(groups), n_loci))
    n_per_group = np.zeros(len(groups))
    for gi, g in enumerate(groups):
        Xg = X[labels == g]
        n_per_group[gi] = Xg.shape[0]
        pg = np.nanmean(Xg, axis=0) / 2.0
        hs_per_group[gi] = 2 * pg * (1 - pg)
    weights = n_per_group / n_per_group.sum()
    hs = np.average(hs_per_group, axis=0, weights=weights)
    return hs, ht, groups


def fst_multi(X, labels):
    hs, ht, groups = group_hs_ht(X, labels)
    ok = np.isfinite(hs) & np.isfinite(ht) & (ht > 0)
    return 1.0 - (np.nanmean(hs[ok]) / np.nanmean(ht[ok])), ok.sum()


def fst_pair(X, labels, ga, gb):
    mask = np.isin(labels, [ga, gb])
    return fst_multi(X[mask], labels[mask])


for panel in PANELS:
    print(f"\n=== [{panel}] Fst (Nei's Gst, Hs/Ht) ===")
    ids, X = load_dosage(panel)

    q_path = TAB / f"admixture_{panel.lower()}_Q.csv"
    q = pd.read_csv(q_path).set_index("sample")
    q = q.reindex(ids)   # align to the dosage matrix's sample order
    qcols = [c for c in q.columns if c.startswith("Q")]
    labels = q[qcols].to_numpy().argmax(axis=1)
    print(f"  n={len(ids)}, {X.shape[1]} SNPs, {len(np.unique(labels))} admixture clusters (elbow K)")

    groups = np.unique(labels)
    pair_rows = []
    for ga, gb in combinations(groups, 2):
        fst, n_loci = fst_pair(X, labels, ga, gb)
        na = int((labels == ga).sum()); nb = int((labels == gb).sum())
        pair_rows.append(dict(panel=panel, cluster_a=int(ga), cluster_b=int(gb),
                               n_a=na, n_b=nb, fst=round(fst, 5), n_loci_used=int(n_loci)))
    pair_df = pd.DataFrame(pair_rows)
    pair_df.to_csv(TAB / f"fst_pairwise_{panel.lower()}.csv", index=False)
    print(f"  pairwise Fst range: [{pair_df['fst'].min():.4f}, {pair_df['fst'].max():.4f}]  "
          f"({len(pair_df)} pairs among {len(groups)} clusters)")

    global_fst, n_loci_g = fst_multi(X, labels)
    perm_fst = np.empty(N_PERM)
    for i in range(N_PERM):
        perm_labels = rng.permutation(labels)
        perm_fst[i], _ = fst_multi(X, perm_labels)
    perm_p = (np.sum(perm_fst >= global_fst) + 1) / (N_PERM + 1)
    print(f"  global multi-cluster Fst = {global_fst:.5f}  (permutation p = {perm_p:.4f}, n_perm={N_PERM})")

    pd.DataFrame([dict(panel=panel, n_clusters=len(groups), n=len(ids),
                        global_fst=round(global_fst, 5), n_loci_used=int(n_loci_g),
                        perm_p=round(perm_p, 4), n_perm=N_PERM)]).to_csv(
        TAB / f"fst_global_{panel.lower()}.csv", index=False)

print("\n-> tables/fst_pairwise_{set}.csv + fst_global_{set}.csv")
