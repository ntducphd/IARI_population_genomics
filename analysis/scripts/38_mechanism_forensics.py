#!/usr/bin/env python
"""38_mechanism_forensics.py — [stage 38] mechanism forensics and robustness checks:

S8  K±1 battery forensics: recompute the k-means K-1/K+1 targets independently, print label
    distributions and confusion vs the elbow labels, and score under BOTH the auto fold rule
    and the primary protocol's fold count — settles whether the byte-identical accuracies in
    S14 were a caching bug or a small-n coincidence, and documents the fold counts.
S1  Size-residualised colour+NIR classifier: residualise each of the 93 colour+NIR features on
    the first 10 PCs of the size/morphology family (OLS, in-sample — standard confound
    removal; residualising on all 111 size features would be n<p and remove everything), then
    re-run the seed-42 classifier protocol. The accuracy drop bounds how much of the
    colour+NIR discriminative signal is stature read through spectral proxies.
S18 Temporal stability: per-imaging-day Mantel (img68/img75/img83 feature blocks) vs genomic
    distance, 999 permutations, both panels.
S15 Pst-Fst bridge: one-way between/within-cluster variance components for family PC1s + the
    height index; Pst(c/h2) sweep 0.1-2.0 compared to the panel's Weir-Cockerham theta.
S7  FDR ledger completion: append the stature-battery, Table-4e partial, and temporal Mantel
    P-values to the Mantel-family BH ledger (robustness_fdr_mantel.csv, feeds Supp Table S12).

Outputs (tables/): kbattery_forensics.csv, residualised_classifier.csv,
temporal_mantel.csv, pst_fst_sweep.csv, robustness_fdr_mantel.csv (extended, overwritten)
"""
import importlib.util
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths import TAB

_spec = importlib.util.spec_from_file_location("stage30", HERE / "30_confound_robustness.py")
s30 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s30)

ELBOW_K = {"Set1": 7, "Set2": 9}


def main():
    kb_rows, rc_rows, tm_rows, pst_rows, fdr_extra = [], [], [], [], []
    for panel in ["Set1", "Set2"]:
        common, coh, Dg, y, _ = s30.load_panel(panel)
        img_cols = [c for c in coh.columns if c.startswith("img")]
        fam = {c: s30.feature_family(c) for c in img_cols}
        size_cols = [c for c in img_cols if fam[c] == "Size_morphology"]
        nonsize_cols = [c for c in img_cols if fam[c] != "Size_morphology"]

        def std(cols):
            X = coh[cols].to_numpy(dtype=float)
            X = np.nan_to_num(X, nan=np.nanmean(X))
            return StandardScaler().fit_transform(X)

        Xs_all, Xs_size, Xs_non = std(img_cols), std(size_cols), std(nonsize_cols)

        # ---- S8: K±1 forensics ----
        pca_t = pd.read_csv(TAB / f"pca_{panel.lower()}.csv")
        pca_t["sample"] = pca_t["sample"].astype(str).str.strip()
        pcs = pca_t.set_index("sample").reindex(common)[
            [c for c in pca_t.columns if c.startswith("PC")]].to_numpy()
        assert np.isfinite(pcs).all(), f"{panel}: NaN in reindexed PCs — ID mismatch"
        yk_store = {}
        for k in (ELBOW_K[panel] - 1, ELBOW_K[panel] + 1):
            yk = KMeans(n_clusters=k, n_init=10, random_state=42).fit_predict(pcs)
            yk_store[k] = yk
            acc_auto, folds_auto, _, _, _ = s30.clf_acc(Xs_all, yk)
            # score under the primary protocol's fold count too (Set1 primary was 2-fold)
            from sklearn.ensemble import RandomForestClassifier
            from sklearn.model_selection import StratifiedKFold, cross_val_predict
            _, counts = np.unique(yk, return_counts=True)
            k2 = max(2, min(2 if panel == "Set1" else 5, counts.min()))
            skf = StratifiedKFold(n_splits=k2, shuffle=True, random_state=42)
            pred2 = cross_val_predict(RandomForestClassifier(
                n_estimators=300, random_state=42, n_jobs=-1), Xs_all, yk, cv=skf)
            kb_rows.append({"panel": panel, "target": f"kmeans_K{k}",
                            "n_correct_auto": int((s30.clf_acc(Xs_all, yk)[4] == yk).sum()),
                            "acc_auto_folds": acc_auto, "folds_auto": folds_auto,
                            "acc_primary_folds": float((pred2 == yk).mean()),
                            "folds_primary": k2,
                            "label_sizes": "/".join(map(str, sorted(counts)))})
            print(f"[38] {panel} kmeans_K{k}: auto {acc_auto:.4f} ({folds_auto}f), "
                  f"primary-protocol {float((pred2 == yk).mean()):.4f} ({k2}f), "
                  f"sizes {sorted(counts)}")
        klo, khi = ELBOW_K[panel] - 1, ELBOW_K[panel] + 1
        same = int((yk_store[klo] == yk_store[khi]).sum())
        print(f"[38] {panel} K-1 vs K+1 labels identical for {same}/{len(common)} accessions "
              f"(different partitions confirmed)" if same < len(common) else
              f"[38] {panel} WARNING: K-1 and K+1 partitions IDENTICAL — bug")

        # ---- S1: size-residualised colour+NIR classifier ----
        p_size = PCA(n_components=10, random_state=42).fit_transform(Xs_size)
        Z = np.column_stack([np.ones(len(p_size)), p_size])
        beta, *_ = np.linalg.lstsq(Z, Xs_non, rcond=None)
        resid = Xs_non - Z @ beta
        resid = StandardScaler().fit_transform(resid)
        acc_r, folds_r, _, _, _ = s30.clf_acc(resid, y)
        acc_o, folds_o, _, _, _ = s30.clf_acc(Xs_non, y)
        rc_rows.append({"panel": panel, "features": "nonsize_residualised_on_sizePC10",
                        "accuracy": acc_r, "cv_folds": folds_r,
                        "nonsize_raw_accuracy": acc_o})
        print(f"[38] {panel} colour+NIR residualised on size-PC10: acc {acc_r:.3f} "
              f"(raw colour+NIR {acc_o:.3f})")

        # ---- S18: per-imaging-day Mantel ----
        for day in ("68", "75", "83"):
            cols = [c for c in img_cols if c.startswith(f"img{day}_")]
            if not cols:
                continue
            Dp_day = s30.dist_std(coh[cols].to_numpy(dtype=float))
            r, p = s30.mantel(Dg, Dp_day)
            tm_rows.append({"panel": panel, "imaging_day_DAS": int(day),
                            "n_features": len(cols), "mantel_r": r, "p": p})
            fdr_extra.append({"panel": panel, "test": f"genomic~phenomic_day{day}", "p": p})
            print(f"[38] {panel} day-{day} Mantel: r = {r:.3f}, P = {p:.3f}")

        # ---- S15: Pst-Fst sweep ----
        theta = pd.read_csv(TAB / f"fst_wc_global_{panel.lower()}.csv").iloc[0]
        fam_pc1 = {}
        for tag, X in (("Size_morphology", Xs_size), ("Colour", std(
                [c for c in img_cols if fam[c] == "Colour"])),
                ("NIR", std([c for c in img_cols if fam[c] == "NIR"]))):
            fam_pc1[tag] = PCA(n_components=1, random_state=42).fit_transform(X)[:, 0]
        height_cols = [c for c in img_cols if ("Height" in c) or re.search(r"_iPH\b", c)]
        H = std(height_cols)
        fam_pc1["height_index"] = H.mean(axis=1)
        for tag, v in fam_pc1.items():
            groups = [v[y == c] for c in np.unique(y) if (y == c).sum() >= 2]
            gm = np.array([g.mean() for g in groups])
            ns = np.array([len(g) for g in groups])
            grand = np.concatenate(groups).mean()
            ssb = float((ns * (gm - grand) ** 2).sum())
            ssw = float(sum(((g - g.mean()) ** 2).sum() for g in groups))
            dfb, dfw = len(groups) - 1, sum(ns) - len(groups)
            msb, msw = ssb / dfb, ssw / dfw
            n0 = (ns.sum() - (ns ** 2).sum() / ns.sum()) / dfb
            s2b = max(0.0, (msb - msw) / n0)
            for ch in (0.1, 0.25, 0.5, 1.0, 2.0):
                pst = ch * s2b / (ch * s2b + 2 * msw) if (s2b + msw) > 0 else 0.0
                pst_rows.append({"panel": panel, "trait": tag, "c_over_h2": ch,
                                 "Pst": pst, "theta_WC": float(theta.theta_global),
                                 "theta_ci_lo": float(theta.ci95_lo),
                                 "theta_ci_hi": float(theta.ci95_hi)})
            p1 = [r_ for r_ in pst_rows if r_["panel"] == panel and r_["trait"] == tag
                  and r_["c_over_h2"] == 1.0][0]["Pst"]
            print(f"[38] {panel} Pst({tag}) at c/h2=1: {p1:.3f} vs theta "
                  f"{theta.theta_global:.3f}")

        # 4e partials + stature battery into the FDR ledger
        cm = pd.read_csv(TAB / f"confound_mantel_{panel.lower()}.csv")
        for _, rr in cm.iterrows():
            fdr_extra.append({"panel": panel, "test": f"confound:{rr['test']}",
                              "p": float(rr["p"])})
        cf = pd.read_csv(TAB / f"concordance_confounder_{panel.lower()}.csv").iloc[0]
        fdr_extra.append({"panel": panel, "test": "phenomic~trait|genomic(partial)",
                          "p": float(cf["p_partial"])})

    pd.DataFrame(kb_rows).to_csv(TAB / "kbattery_forensics.csv", index=False)
    pd.DataFrame(rc_rows).to_csv(TAB / "residualised_classifier.csv", index=False)
    pd.DataFrame(tm_rows).to_csv(TAB / "temporal_mantel.csv", index=False)
    pd.DataFrame(pst_rows).to_csv(TAB / "pst_fst_sweep.csv", index=False)

    # ---- S7: extend the Mantel-family BH ledger ----
    led = pd.read_csv(TAB / "robustness_fdr_mantel.csv")
    keep_cols = list(led.columns)
    add = pd.DataFrame(fdr_extra)
    add = add[~add.apply(lambda r_: ((led.get("test") == r_["test"]) &
                                     (led.get("panel") == r_["panel"])).any(), axis=1)]
    merged = pd.concat([led[["panel", "test", "p"]] if set(
        ["panel", "test", "p"]).issubset(led.columns) else led, add], ignore_index=True)
    m = merged.sort_values("p").reset_index(drop=True)
    n = len(m)
    q = (m["p"] * n / (np.arange(n) + 1))[::-1].cummin()[::-1]
    m["bh_q"] = q
    m["significant_q05"] = m["bh_q"] < 0.05
    m.to_csv(TAB / "robustness_fdr_mantel.csv", index=False)
    print(f"[38] FDR ledger extended: {len(led)} -> {n} tests, "
          f"{int(m['significant_q05'].sum())} survive q<0.05")
    print("[38] done")


if __name__ == "__main__":
    main()
