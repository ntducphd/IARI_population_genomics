#!/usr/bin/env python
"""37_mechanism_stature_tests.py — [stage 37] mechanism tests for the concordance signal, computable from data on
disk (review_6_nature_round2.md, R2 major comments):

(a) SIZE-LEAKAGE TEST for the colour+NIR subset. The manuscript claimed the colour+NIR
    features "contain no size information" — an assertion, never a measurement. Measured here
    as the cross-validated R^2 (5-fold, out-of-sample) of predicting the height index (mean
    z-score of the 9 height features, stage-36 definition) from the 93 colour+NIR features,
    with both a linear (ridge, alpha=1) and a random-forest regressor. High CV-R^2 would mean
    the "more than a ruler" claim needs bounding; the number goes in the paper either way
    (không bịa).

(b) CLASSIFIER UNCERTAINTY. All subset-accuracy comparisons were point estimates on a single
    CV split. Here: 20 repeats of the identical stratified-CV protocol with different fold
    seeds for every (panel x feature-subset) cell and the external-label test -> mean, SD,
    and 2.5/97.5 percentiles across repeats. The original seed-42 point estimate stays the
    headline (protocol unchanged); the spread is reported alongside it.

Outputs (analysis/results/tables/): sizeleakage_summary.csv, classifier_uncertainty.csv
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold, StratifiedKFold, cross_val_predict, cross_val_score
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths import TAB

_spec = importlib.util.spec_from_file_location("stage30", HERE / "30_confound_robustness.py")
s30 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s30)

import re

N_REPEATS = 20


def repeat_cv_acc(X, y, n_repeats=N_REPEATS):
    _, counts = np.unique(y, return_counts=True)
    k = max(2, min(5, counts.min()))
    accs = []
    for rep in range(n_repeats):
        skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=2000 + rep)
        clf = RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=-1)
        pred = cross_val_predict(clf, X, y, cv=skf)
        accs.append(float((pred == y).mean()))
    a = np.array(accs)
    return k, float(a.mean()), float(a.std(ddof=1)), float(np.percentile(a, 2.5)), \
        float(np.percentile(a, 97.5))


def main():
    leak_rows, unc_rows = [], []
    for panel in ["Set1", "Set2"]:
        common, coh, Dg, y, _ = s30.load_panel(panel)
        img_cols = [c for c in coh.columns if c.startswith("img")]
        fam = {c: s30.feature_family(c) for c in img_cols}
        nonsize_cols = [c for c in img_cols if fam[c] != "Size_morphology"]
        height_cols = [c for c in img_cols
                       if ("Height" in c) or re.search(r"_iPH\b", c)]

        def std(cols):
            X = coh[cols].to_numpy(dtype=float)
            X = np.nan_to_num(X, nan=np.nanmean(X))
            return StandardScaler().fit_transform(X)

        Xs_all, Xs_ht, Xs_non = std(img_cols), std(height_cols), std(nonsize_cols)
        H = coh[height_cols].to_numpy(dtype=float)
        H = np.nan_to_num(H, nan=np.nanmean(H))
        Hz = (H - H.mean(axis=0)) / H.std(axis=0)
        hidx = Hz.mean(axis=1)

        # ---- (a) size leakage: predict height index from colour+NIR, out of sample ----
        cv = KFold(n_splits=5, shuffle=True, random_state=42)
        for name, model in [("ridge", Ridge(alpha=1.0)),
                            ("random_forest", RandomForestRegressor(
                                n_estimators=300, random_state=42, n_jobs=-1))]:
            r2 = cross_val_score(model, Xs_non, hidx, cv=cv, scoring="r2")
            leak_rows.append({"panel": panel, "model": name,
                              "cv_r2_mean": float(r2.mean()),
                              "cv_r2_sd": float(r2.std(ddof=1)),
                              "n_features": len(nonsize_cols), "n": len(hidx)})
            print(f"[37] {panel} height~colour+NIR leakage ({name}): "
                  f"CV-R2 = {r2.mean():.3f} ± {r2.std(ddof=1):.3f}")

        # ---- (b) repeated-CV uncertainty for the classifier battery ----
        for tag, Xf in [("full_204", Xs_all), ("height_only", Xs_ht),
                        ("nonsize_colour_nir", Xs_non)]:
            k, m, sd, lo, hi = repeat_cv_acc(Xf, y)
            unc_rows.append({"panel": panel, "features": tag,
                             "target": "admixture_elbowK", "cv_folds": k,
                             "n_repeats": N_REPEATS, "acc_mean": m, "acc_sd": sd,
                             "acc_p2.5": lo, "acc_p97.5": hi})
            print(f"[37] {panel} clf {tag}: {m:.3f} ± {sd:.3f} [{lo:.3f}, {hi:.3f}]")

        if panel == "Set1":
            sub = pd.read_csv(TAB / "subpop_assignment_set1.csv")
            m_ = dict(zip(sub["sample"].astype(str), sub["subpopulation"].astype(str)))
            lab = pd.Series([m_.get(s, "NA") for s in common])
            counts = lab.value_counts()
            keepg = lab.isin([g for g in counts.index if g != "NA" and counts[g] >= 5])
            ye, _names = pd.factorize(lab[keepg])
            Xe = Xs_all[keepg.to_numpy()]
            k, m, sd, lo, hi = repeat_cv_acc(Xe, ye)
            unc_rows.append({"panel": "Set1", "features": "full_204",
                             "target": "external_3KRGP_subpop", "cv_folds": k,
                             "n_repeats": N_REPEATS, "acc_mean": m, "acc_sd": sd,
                             "acc_p2.5": lo, "acc_p97.5": hi})
            print(f"[37] Set1 clf external: {m:.3f} ± {sd:.3f} [{lo:.3f}, {hi:.3f}]")

    pd.DataFrame(leak_rows).to_csv(TAB / "sizeleakage_summary.csv", index=False)
    pd.DataFrame(unc_rows).to_csv(TAB / "classifier_uncertainty.csv", index=False)
    print("[37] done")


if __name__ == "__main__":
    main()
