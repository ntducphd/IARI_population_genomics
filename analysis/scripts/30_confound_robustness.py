#!/usr/bin/env python
"""30_confound_robustness.py — [stage 30] confound and label-robustness batteries for the concordance analysis
explicit confound tests and label-robustness for the Pillar B claim.

[F1] Size/stature confounding:
  * D_height  = distance on the height features only (HeightSV/HeightTV/iPH x 3 timepoints)
  * D_size    = distance on the full size/morphology family (111 features)
  * D_nonsize = distance on colour + NIR features only (93 features)
  * partial Mantel r(Dg, Dp | D_size)  — does ANY full-set signal survive controlling stature?
  * Mantel r(Dg, D_height) and r(Dg, D_nonsize | D_size) — how far does height alone go, and is
    there residual non-size signal?
  * Height-only classifier baseline vs full-204 (same CV protocol) — the deflationary test.
[F1-spatial] PlantID is a numeric pot/sequence identifier; |ΔPlantID| is used as a greenhouse-
  position PROXY: Mantel r(Dg, D_pos), r(Dp, D_pos), and partial r(Dg, Dp | D_pos). (An honest
  proxy, stated as such; true coordinates were not recorded in the data lake.)
[F3] Bootstrap 95% CIs (500 resamples over accessions) for the primary Mantel r, both panels;
  primary classifier empirical P upgraded to 999 permutations.
[F4] Label-robustness for the classifier: (a) EXTERNAL 3K-RGP subpopulation labels as target
  (Set 1, groups n >= 5) — breaks the circularity objection; (b) k-means labels at K-1/K/K+1;
  (c) soft-Q-weighted accuracy (weight = accession's max admixture proportion).

Outputs (analysis/results/tables/): confound_mantel_{set}.csv, classifier_baselines_{set}.csv,
external_label_classifier_set1.csv, spatial_proxy_{set}.csv, mantel_bootstrap_{set}.csv,
stage30_summary.csv
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import TAB, WORKSPACE

COHORT = {
    "Set1": WORKSPACE / "manuscript_7_phenomic_selection/analysis/data/input/cohort_set1.csv",
    "Set2": WORKSPACE / "manuscript_7_phenomic_selection/analysis/data/input/cohort_set2.csv",
}
ELBOW_K = {"Set1": 7, "Set2": 9}
N_PERM = 999
N_PERM_CLF = 999
N_BOOT = 500
RNG = np.random.default_rng(42)


def normalize_id(x):
    return re.sub(r"[\s\-]", "", str(x).upper())


def feature_family(col):
    if "Color." in col:
        return "Colour"
    if col.startswith(("img68_NIR", "img75_NIR", "img83_NIR",
                       "img68_IR.", "img75_IR.", "img83_IR.")):
        return "NIR"
    return "Size_morphology"


def upper(d):
    return d[np.triu_indices_from(d, k=1)]


def mantel(d1, d2, n_perm=N_PERM):
    u1, u2 = upper(d1), upper(d2)
    r_obs, _ = pearsonr(u1, u2)
    count = 0
    n = d1.shape[0]
    for _ in range(n_perm):
        p = RNG.permutation(n)
        r_p, _ = pearsonr(u1, upper(d2[np.ix_(p, p)]))
        count += r_p >= r_obs
    return float(r_obs), (count + 1) / (n_perm + 1)


def partial_mantel(d1, d2, d3, n_perm=N_PERM):
    u1, u2, u3 = upper(d1), upper(d2), upper(d3)

    def pr(a, b, c):
        rab, _ = pearsonr(a, b)
        rac, _ = pearsonr(a, c)
        rbc, _ = pearsonr(b, c)
        return (rab - rac * rbc) / np.sqrt((1 - rac**2) * (1 - rbc**2))

    r_obs = pr(u1, u2, u3)
    count = 0
    n = d1.shape[0]
    for _ in range(n_perm):
        p = RNG.permutation(n)
        count += pr(u1, upper(d2[np.ix_(p, p)]), u3) >= r_obs
    return float(r_obs), (count + 1) / (n_perm + 1)


def dist_std(X):
    X = np.nan_to_num(X, nan=np.nanmean(X))
    return squareform(pdist(StandardScaler().fit_transform(X), metric="euclidean"))


def cv_k(y):
    _, counts = np.unique(y, return_counts=True)
    return max(2, min(5, counts.min()))


def clf_acc(X, y, n_perm=0, rng=RNG):
    k = cv_k(y)
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=42)
    clf = RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=-1)
    pred = cross_val_predict(clf, X, y, cv=skf)
    acc = float((pred == y).mean())
    p_emp, null_mean = None, None
    if n_perm:
        null = []
        for i in range(n_perm):
            yp = rng.permutation(y)
            try:
                null.append((cross_val_predict(
                    clf, X, yp, cv=StratifiedKFold(n_splits=k, shuffle=True,
                                                   random_state=1000 + i)) == yp).mean())
            except ValueError:
                continue
        null = np.array(null)
        p_emp = float((np.sum(null >= acc) + 1) / (len(null) + 1))
        null_mean = float(null.mean())
    return acc, k, p_emp, null_mean, pred


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
    q = pd.read_csv(TAB / f"admixture_{panel.lower()}_Q.csv")
    q["sample"] = q["sample"].astype(str).str.strip()
    qcols = [c for c in q.columns if c.startswith("Q")]
    qmap = dict(zip(q["sample"], q[qcols].to_numpy().argmax(axis=1)))
    softq = dict(zip(q["sample"], q[qcols].to_numpy().max(axis=1)))
    y = np.array([qmap.get(s, -1) for s in common])
    w = np.array([softq.get(s, np.nan) for s in common])
    keep = y >= 0
    return (np.asarray(common)[keep], coh.loc[np.asarray(common)[keep]],
            Dg[np.ix_(keep, keep)], y[keep], w[keep])


def main():
    summary = []
    for panel in ["Set1", "Set2"]:
        common, coh, Dg, y, softw = load_panel(panel)
        img_cols = [c for c in coh.columns if c.startswith("img")]
        fam = {c: feature_family(c) for c in img_cols}
        size_cols = [c for c in img_cols if fam[c] == "Size_morphology"]
        nonsize_cols = [c for c in img_cols if fam[c] != "Size_morphology"]
        height_cols = [c for c in img_cols
                       if ("Height" in c) or re.search(r"_iPH\b", c)]
        X_all = coh[img_cols].to_numpy(dtype=float)
        Dp = dist_std(X_all)
        Dsize = dist_std(coh[size_cols].to_numpy(dtype=float))
        Dnon = dist_std(coh[nonsize_cols].to_numpy(dtype=float))
        Dht = dist_std(coh[height_cols].to_numpy(dtype=float))
        print(f"[30] {panel}: n={len(common)}, height features = {len(height_cols)}")

        # ---- F1 Mantel battery ----
        rows = []
        for name, args, partial in [
            ("genomic~phenomic (primary)", (Dg, Dp), None),
            ("genomic~height_only", (Dg, Dht), None),
            ("genomic~phenomic | size_family", (Dg, Dp, Dsize), "partial"),
            ("genomic~nonsize | size_family", (Dg, Dnon, Dsize), "partial"),
            ("genomic~size_family", (Dg, Dsize), None),
        ]:
            if partial:
                r, p = partial_mantel(*args)
            else:
                r, p = mantel(*args)
            rows.append({"panel": panel, "test": name, "r": r, "p": p})
            print(f"    {name}: r = {r:.3f}, P = {p:.3f}")
        pd.DataFrame(rows).to_csv(TAB / f"confound_mantel_{panel.lower()}.csv", index=False)

        # ---- F1 spatial proxy (PlantID sequence) ----
        pid = pd.to_numeric(coh["PlantID"], errors="coerce").to_numpy(dtype=float)
        sp_rows = []
        if np.isfinite(pid).sum() > 0.9 * len(pid):
            Dpos = np.abs(pid[:, None] - pid[None, :])
            for name, args, partial in [
                ("genomic~position_proxy", (Dg, Dpos), None),
                ("phenomic~position_proxy", (Dp, Dpos), None),
                ("genomic~phenomic | position_proxy", (Dg, Dp, Dpos), "partial"),
            ]:
                if partial:
                    r, p = partial_mantel(*args)
                else:
                    r, p = mantel(*args)
                sp_rows.append({"panel": panel, "test": name, "r": r, "p": p})
                print(f"    {name}: r = {r:.3f}, P = {p:.3f}")
        else:
            sp_rows.append({"panel": panel, "test": "PlantID not numeric", "r": np.nan,
                            "p": np.nan})
        pd.DataFrame(sp_rows).to_csv(TAB / f"spatial_proxy_{panel.lower()}.csv", index=False)

        # ---- F3 bootstrap CI for the primary Mantel r ----
        n = Dg.shape[0]
        boots = []
        for _ in range(N_BOOT):
            idx = RNG.integers(0, n, n)
            idx = np.unique(idx)                      # unique to keep distances meaningful
            if len(idx) < 20:
                continue
            r_b, _ = pearsonr(upper(Dg[np.ix_(idx, idx)]), upper(Dp[np.ix_(idx, idx)]))
            boots.append(r_b)
        lo, hi = np.percentile(boots, [2.5, 97.5])
        pd.DataFrame([{"panel": panel, "r_boot_mean": float(np.mean(boots)),
                       "ci95_lo": float(lo), "ci95_hi": float(hi),
                       "n_boot": len(boots)}]).to_csv(
            TAB / f"mantel_bootstrap_{panel.lower()}.csv", index=False)
        print(f"    bootstrap CI (primary r): [{lo:.3f}, {hi:.3f}]")

        # ---- F1/F4 classifier battery ----
        Xs_all = StandardScaler().fit_transform(np.nan_to_num(X_all, nan=np.nanmean(X_all)))
        Xs_ht = StandardScaler().fit_transform(
            np.nan_to_num(coh[height_cols].to_numpy(dtype=float),
                          nan=np.nanmean(coh[height_cols].to_numpy(dtype=float))))
        Xs_non = StandardScaler().fit_transform(
            np.nan_to_num(coh[nonsize_cols].to_numpy(dtype=float),
                          nan=np.nanmean(coh[nonsize_cols].to_numpy(dtype=float))))
        crows = []
        acc_full, k_used, p_full, null_full, pred_full = clf_acc(Xs_all, y,
                                                                 n_perm=N_PERM_CLF)
        soft_acc = float(np.average(pred_full == y, weights=softw))
        crows.append({"panel": panel, "features": "full_204", "target": "admixture_elbowK",
                      "accuracy": acc_full, "cv_folds": k_used, "empirical_p": p_full,
                      "null_mean": null_full, "softQ_weighted_accuracy": soft_acc})
        for tag, Xf in [("height_only", Xs_ht), ("nonsize_colour_nir", Xs_non)]:
            a, k2, _, _, _ = clf_acc(Xf, y)
            crows.append({"panel": panel, "features": tag, "target": "admixture_elbowK",
                          "accuracy": a, "cv_folds": k2, "empirical_p": None,
                          "null_mean": None, "softQ_weighted_accuracy": None})
        # K±1 k-means labels (PCA coords)
        pca = pd.read_csv(TAB / f"pca_{panel.lower()}.csv")
        pca["sample"] = pca["sample"].astype(str).str.strip()
        pcs = pca.set_index("sample").reindex(common)[
            [c for c in pca.columns if c.startswith("PC")]].to_numpy()
        for k in (ELBOW_K[panel] - 1, ELBOW_K[panel] + 1):
            yk = KMeans(n_clusters=k, n_init=10, random_state=42).fit_predict(pcs)
            a, k2, _, _, _ = clf_acc(Xs_all, yk)
            crows.append({"panel": panel, "features": "full_204",
                          "target": f"kmeans_K{k}", "accuracy": a, "cv_folds": k2,
                          "empirical_p": None, "null_mean": None,
                          "softQ_weighted_accuracy": None})
        pd.DataFrame(crows).to_csv(TAB / f"classifier_baselines_{panel.lower()}.csv",
                                   index=False)
        for r_ in crows:
            print(f"    clf [{r_['features']} -> {r_['target']}]: acc = {r_['accuracy']:.3f}")

        # ---- F4 external labels (Set 1 only) ----
        if panel == "Set1":
            sub = pd.read_csv(TAB / "subpop_assignment_set1.csv")
            m = dict(zip(sub["sample"].astype(str), sub["subpopulation"].astype(str)))
            lab = pd.Series([m.get(s, "NA") for s in common])
            counts = lab.value_counts()
            keepg = lab.isin([g for g in counts.index if g != "NA" and counts[g] >= 5])
            ye, Xe = pd.factorize(lab[keepg])[0], Xs_all[keepg.to_numpy()]
            acc_e, ke, pe, nulle, _ = clf_acc(Xe, ye, n_perm=N_PERM_CLF)
            pd.DataFrame([{"panel": "Set1", "target": "external_3KRGP_subpop",
                           "n": int(keepg.sum()), "n_groups": int(lab[keepg].nunique()),
                           "accuracy": acc_e, "cv_folds": ke, "empirical_p": pe,
                           "null_mean": nulle}]).to_csv(
                TAB / "external_label_classifier_set1.csv", index=False)
            print(f"    clf [full_204 -> EXTERNAL 3K-RGP]: acc = {acc_e:.3f} "
                  f"(null {nulle:.3f}, P = {pe:.4f})")

        summary.append({"panel": panel, "n": len(common),
                        "n_height_features": len(height_cols)})
    pd.DataFrame(summary).to_csv(TAB / "stage30_summary.csv", index=False)
    print("[30] done")


if __name__ == "__main__":
    main()
