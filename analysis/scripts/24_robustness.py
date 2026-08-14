#!/usr/bin/env python
"""24_robustness.py — [stage 24] robustness layer for the Pillar B concordance claims.

Motivation: the core inference
leans on the Mantel family, which is the field's most criticised test (Guillot & Rousset 2013).
This stage adds the defence lines a methods referee will ask for:

  [A] Classifier nulls: same protocol as stage 09 (random forest, 300 trees, stratified
      k-fold with k = largest value <= 5 giving every cluster >= k members), now reporting
      balanced accuracy, macro-F1, per-class recall, and an EMPIRICAL label-permutation null
      (N_PERM_CLF full CV repeats with shuffled cluster labels -> empirical P for the observed
      accuracy), replacing the majority-class baseline as the formal chance reference.
  [B] MMRR: multiple matrix regression with randomisation (Wang 2013) — phenomic distance
      ~ genomic distance + trait distance, coefficients tested by row/column permutation of the
      response matrix. A regression-on-distance-matrices alternative to (partial) Mantel.
  [C] Procrustes PC-count sensitivity: M^2 and PROTEST P for n_pc = 2..8 (stage 09 fixed
      n_pc = 4).
  [D] Elbow-K sensitivity: ARI between PCA k-means at K-1/K/K+1 and the admixture argmax
      assignment (elbow K fixed) — does the cross-method consensus survive +/-1 in K?
  [E] FDR ledger: every Mantel-family P reported in the paper, Benjamini-Hochberg adjusted
      in one table.

Outputs (analysis/results/tables/): robustness_classifier_{set}.csv,
robustness_mmrr_{set}.csv, robustness_procrustes_{set}.csv, robustness_ari_{set}.csv,
robustness_fdr_mantel.csv, robustness_summary.csv
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import procrustes as scipy_procrustes
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (adjusted_rand_score, balanced_accuracy_score, f1_score,
                             recall_score)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import TAB, COHORT_SET1, COHORT_SET2

COHORT = {"Set1": COHORT_SET1, "Set2": COHORT_SET2}
ELBOW_K = {"Set1": 7, "Set2": 9}
N_PERM_CLF = 100
N_PERM_MMRR = 999
RNG = np.random.default_rng(42)


def normalize_id(x):
    return re.sub(r"[\s\-]", "", str(x).upper())


def load_panel(panel):
    cohort = pd.read_csv(COHORT[panel])
    cohort["Taxa"] = cohort["Taxa"].astype(str).str.strip()
    cohort["_key"] = cohort["Taxa"].map(normalize_id)

    ibs = pd.read_csv(TAB / f"ibs_dist_{panel.lower()}.csv")
    genomic_ids = sorted(set(ibs["sample_a"]) | set(ibs["sample_b"]))
    g_key = {normalize_id(g): g for g in genomic_ids}
    common_keys = sorted(set(cohort["_key"]) & set(g_key))
    common = [g_key[k] for k in common_keys]
    coh = cohort.drop_duplicates("_key").set_index("_key").loc[common_keys]
    coh.index = common

    ibs_p = ibs.pivot(index="sample_a", columns="sample_b", values="distance")
    Dg = ibs_p.reindex(index=common, columns=common).to_numpy().copy()
    np.fill_diagonal(Dg, 0.0)

    img_cols = [c for c in coh.columns if c.startswith("img")]
    X = coh[img_cols].to_numpy(dtype=float)
    X = np.nan_to_num(X, nan=np.nanmean(X))
    Xs = StandardScaler().fit_transform(X)
    Dp = squareform(pdist(Xs, metric="euclidean"))

    trait_cols = [c for c in coh.columns if c.endswith("_Control")]
    Xt = coh[trait_cols].to_numpy(dtype=float)
    cm = np.nanmean(Xt, axis=0)
    idx = np.where(np.isnan(Xt))
    Xt[idx] = np.take(cm, idx[1])
    Dt = squareform(pdist(StandardScaler().fit_transform(Xt), metric="euclidean"))

    q = pd.read_csv(TAB / f"admixture_{panel.lower()}_Q.csv")
    q["sample"] = q["sample"].astype(str).str.strip()
    qcols = [c for c in q.columns if c.startswith("Q")]
    qmap = dict(zip(q["sample"], q[qcols].to_numpy().argmax(axis=1)))
    y = np.array([qmap.get(s, -1) for s in common])
    keep = y >= 0
    return (np.asarray(common)[keep], Xs[keep], Dg[np.ix_(keep, keep)],
            Dp[np.ix_(keep, keep)], Dt[np.ix_(keep, keep)], y[keep])


def cv_k(y):
    _, counts = np.unique(y, return_counts=True)
    return max(2, min(5, counts.min()))


def classifier_block(panel, Xs, y):
    k = cv_k(y)
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=42)
    clf = RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=-1)
    pred = cross_val_predict(clf, Xs, y, cv=skf)
    acc = float((pred == y).mean())
    bal = float(balanced_accuracy_score(y, pred))
    mf1 = float(f1_score(y, pred, average="macro"))
    rec = recall_score(y, pred, average=None, labels=np.unique(y))

    null = []
    for i in range(N_PERM_CLF):
        yp = RNG.permutation(y)
        try:
            null.append((cross_val_predict(clf, Xs, yp,
                        cv=StratifiedKFold(n_splits=k, shuffle=True,
                                           random_state=100 + i)) == yp).mean())
        except ValueError:
            continue
    null = np.array(null)
    p_emp = float((np.sum(null >= acc) + 1) / (len(null) + 1))

    rows = [{"metric": "accuracy", "value": acc},
            {"metric": "balanced_accuracy", "value": bal},
            {"metric": "macro_f1", "value": mf1},
            {"metric": "null_mean_accuracy", "value": float(null.mean())},
            {"metric": "null_p97_5", "value": float(np.percentile(null, 97.5))},
            {"metric": "empirical_p", "value": p_emp},
            {"metric": "cv_folds", "value": k}]
    rows += [{"metric": f"recall_C{c}", "value": float(r)}
             for c, r in zip(np.unique(y), rec)]
    pd.DataFrame(rows).to_csv(TAB / f"robustness_classifier_{panel.lower()}.csv", index=False)
    print(f"[24A] {panel}: acc={acc:.3f} (null {null.mean():.3f}, P={p_emp:.3f}), "
          f"balanced={bal:.3f}, macroF1={mf1:.3f}")
    return acc, p_emp


def mmrr(panel, Dp, Dg, Dt):
    iu = np.triu_indices_from(Dp, k=1)
    Y = (Dp[iu] - Dp[iu].mean()) / Dp[iu].std()
    X = np.column_stack([(D[iu] - D[iu].mean()) / D[iu].std() for D in (Dg, Dt)])
    Xd = np.column_stack([np.ones(len(Y)), X])
    beta = np.linalg.lstsq(Xd, Y, rcond=None)[0]
    n = Dp.shape[0]

    def betas_for(Dperm):
        Yp = (Dperm[iu] - Dperm[iu].mean()) / Dperm[iu].std()
        return np.linalg.lstsq(Xd, Yp, rcond=None)[0]

    count = np.zeros(3)
    for _ in range(N_PERM_MMRR):
        perm = RNG.permutation(n)
        bp = betas_for(Dp[np.ix_(perm, perm)])
        count += (np.abs(bp) >= np.abs(beta)).astype(int)
    pvals = (count + 1) / (N_PERM_MMRR + 1)
    out = pd.DataFrame({"term": ["intercept", "genomic", "trait"],
                        "beta_std": beta, "p_perm": pvals})
    out.to_csv(TAB / f"robustness_mmrr_{panel.lower()}.csv", index=False)
    print(f"[24B] {panel}: MMRR beta_genomic={beta[1]:.3f} (P={pvals[1]:.3f}), "
          f"beta_trait={beta[2]:.3f} (P={pvals[2]:.3f})")


def procrustes_sens(panel, Xs, Dg):
    # classical MDS of the genomic distance for genomic coordinates
    n = Dg.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (Dg ** 2) @ J
    evals, evecs = np.linalg.eigh(B)
    order = np.argsort(evals)[::-1]
    Gfull = evecs[:, order] * np.sqrt(np.maximum(evals[order], 0))
    from sklearn.decomposition import PCA
    Pfull = PCA(n_components=8, random_state=42).fit_transform(Xs)
    rows = []
    for npc in range(2, 9):
        m1, m2, disp = scipy_procrustes(Gfull[:, :npc], Pfull[:, :npc])
        count = 0
        for _ in range(199):
            perm = RNG.permutation(n)
            _, _, dp_ = scipy_procrustes(Gfull[:, :npc], Pfull[perm, :npc])
            count += dp_ <= disp
        rows.append({"n_pc": npc, "m2": float(disp), "p_perm": (count + 1) / 200})
    pd.DataFrame(rows).to_csv(TAB / f"robustness_procrustes_{panel.lower()}.csv", index=False)
    m2s = [f"{r['m2']:.3f}" for r in rows]
    print(f"[24C] {panel}: M2 across n_pc 2-8 = {', '.join(m2s)}")


def ari_sens(panel):
    pca = pd.read_csv(TAB / f"pca_{panel.lower()}.csv")
    q = pd.read_csv(TAB / f"admixture_{panel.lower()}_Q.csv")
    qcols = [c for c in q.columns if c.startswith("Q")]
    merged = pca.merge(q, on="sample")
    coords = merged[[c for c in pca.columns if c.startswith("PC")]].to_numpy()
    admix = merged[qcols].to_numpy().argmax(axis=1)
    K = ELBOW_K[panel]
    rows = []
    for k in (K - 1, K, K + 1):
        km = KMeans(n_clusters=k, n_init=10, random_state=42).fit_predict(coords)
        rows.append({"kmeans_k": k, "ari_vs_admixture_elbowK": float(
            adjusted_rand_score(admix, km))})
    pd.DataFrame(rows).to_csv(TAB / f"robustness_ari_{panel.lower()}.csv", index=False)
    print(f"[24D] {panel}: ARI at K-1/K/K+1 = "
          + ", ".join(f"{r['ari_vs_admixture_elbowK']:.3f}" for r in rows))


def fdr_ledger():
    rows = []
    for panel in ["Set1", "Set2"]:
        m = pd.read_csv(TAB / f"concordance_mantel_{panel.lower()}.csv")
        for _, r in m.iterrows():
            rows.append({"panel": panel, "test": r["comparison"], "p": r["p"]})
        fa = pd.read_csv(TAB / f"concordance_feature_attribution_{panel.lower()}.csv")
        pcol = [c for c in fa.columns if c.lower() in ("p", "pvalue", "p_value")][0]
        gcol = [c for c in fa.columns if "feat" in c.lower() or "group" in c.lower()][0]
        for _, r in fa.iterrows():
            rows.append({"panel": panel, "test": f"family:{r[gcol]}", "p": r[pcol]})
    df = pd.DataFrame(rows)
    order = df["p"].rank(method="first")
    m = len(df)
    df["q_bh"] = np.minimum.accumulate(
        (df["p"] * m / order).iloc[np.argsort(-order.values)].values)[np.argsort(
            np.argsort(-order.values))]
    df["q_bh"] = df["q_bh"].clip(upper=1.0)
    df["significant_q05"] = df["q_bh"] < 0.05
    df.to_csv(TAB / "robustness_fdr_mantel.csv", index=False)
    print(f"[24E] FDR ledger: {df['significant_q05'].sum()}/{len(df)} tests survive BH q<0.05")


def main():
    summary = []
    for panel in ["Set1", "Set2"]:
        common, Xs, Dg, Dp, Dt, y = load_panel(panel)
        print(f"[24] {panel}: n={len(common)} matched accessions")
        acc, p_emp = classifier_block(panel, Xs, y)
        mmrr(panel, Dp, Dg, Dt)
        procrustes_sens(panel, Xs, Dg)
        ari_sens(panel)
        summary.append({"panel": panel, "n": len(common),
                        "clf_accuracy": acc, "clf_empirical_p": p_emp})
    fdr_ledger()
    pd.DataFrame(summary).to_csv(TAB / "robustness_summary.csv", index=False)
    print("[24] done -> robustness_summary.csv")


if __name__ == "__main__":
    main()
