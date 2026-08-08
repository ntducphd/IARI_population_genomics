#!/usr/bin/env python
"""38b_residualised_perm.py — [stage 38b] empirical permutation P for the size-residualised
colour+NIR classifier (stage 38, S1). The residualised accuracy (0.260 Set1 / 0.238 Set2) is
now a load-bearing number for the weakened "size-independent residual" claim, so it gets the
same 999-label-permutation null as every other load-bearing classifier accuracy.
Output: residualised_classifier.csv (overwritten with empirical_p / null_mean columns added).
"""
import importlib.util
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths import TAB

_spec = importlib.util.spec_from_file_location("stage30", HERE / "30_confound_robustness.py")
s30 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s30)

rows = []
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

    Xs_size, Xs_non = std(size_cols), std(nonsize_cols)
    p_size = PCA(n_components=10, random_state=42).fit_transform(Xs_size)
    Z = np.column_stack([np.ones(len(p_size)), p_size])
    beta, *_ = np.linalg.lstsq(Z, Xs_non, rcond=None)
    resid = StandardScaler().fit_transform(Xs_non - Z @ beta)
    acc, folds, p_emp, null_mean, _ = s30.clf_acc(resid, y, n_perm=999)
    rows.append({"panel": panel, "features": "nonsize_residualised_on_sizePC10",
                 "accuracy": acc, "cv_folds": folds, "empirical_p": p_emp,
                 "null_mean": null_mean,
                 "nonsize_raw_accuracy": float(pd.read_csv(
                     TAB / "residualised_classifier.csv").set_index("panel")
                     .loc[panel, "nonsize_raw_accuracy"])})
    print(f"[38b] {panel}: residualised acc {acc:.3f} ({folds}f), null {null_mean:.3f}, "
          f"P = {p_emp:.3f}")

pd.DataFrame(rows).to_csv(TAB / "residualised_classifier.csv", index=False)
print("[38b] done")
