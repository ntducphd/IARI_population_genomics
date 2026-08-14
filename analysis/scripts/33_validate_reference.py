#!/usr/bin/env python
"""33_validate_reference.py — [stage 33] cross-validation of bespoke implementations:
(a) validate the bespoke Mantel implementation against the reference implementation in R vegan
    (same distance matrices exported and tested with vegan::mantel, 999 permutations);
(b) RECONCILE the Procrustes pipelines: the primary analysis (stage 09) used the SNPRelate
    genomic PCs (tables/pca_{set}.csv) vs the sklearn PCA of standardised imaging features,
    while the stage-24 sensitivity used a classical-MDS reconstruction of the IBS distance —
    a different genomic configuration, which produced the M^2 clash the reviews flagged
    (0.805/0.856 primary vs 0.864/0.888 at the same 4 PCs). This stage recomputes the
    sensitivity across n_pc = 2..6 with EXACTLY the stage-09 inputs, so the primary value sits
    inside its own band, and overwrites robustness_procrustes_{set}.csv.
    (pca_{set}.csv provides 6 PCs, so the band is 2-6, stated in the manuscript.)

Outputs: mantel_vegan_validation.csv, robustness_procrustes_{set}.csv (recomputed),
interim/valid_{set}_{dg,dp}.csv (exported matrices for the R step)
"""
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import procrustes as scipy_procrustes
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import TAB, INTERIM, RSCRIPT, COHORT_SET1, COHORT_SET2

COHORT = {"Set1": COHORT_SET1, "Set2": COHORT_SET2}
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
    return common, Xs, Dg, Dp


R_VALIDATE = r'''
suppressMessages(library(vegan))
args <- commandArgs(trailingOnly = TRUE)
dg <- as.dist(as.matrix(read.csv(args[1], row.names = 1)))
dp <- as.dist(as.matrix(read.csv(args[2], row.names = 1)))
set.seed(42)
m <- vegan::mantel(dg, dp, permutations = 999, method = "pearson")
cat(sprintf("VEGAN %s r=%.6f p=%.4f\n", args[3], m$statistic, m$signif))
'''


def main():
    rows = []
    for panel in ["Set1", "Set2"]:
        common, Xs, Dg, Dp = load_panel(panel)

        # ---- (a) export + vegan Mantel ----
        dg_path = INTERIM / f"valid_{panel.lower()}_dg.csv"
        dp_path = INTERIM / f"valid_{panel.lower()}_dp.csv"
        pd.DataFrame(Dg, index=common, columns=common).to_csv(dg_path)
        pd.DataFrame(Dp, index=common, columns=common).to_csv(dp_path)
        rscript = INTERIM / "validate_mantel.R"
        rscript.write_text(R_VALIDATE)
        r = subprocess.run([str(RSCRIPT), str(rscript), str(dg_path), str(dp_path), panel],
                           capture_output=True, text=True)
        out = (r.stdout or "") + (r.stderr or "")
        mm = re.search(r"VEGAN \S+ r=([\d.\-]+) p=([\d.]+)", out)
        if not mm:
            print(f"[33a] {panel}: vegan run failed:\n{out[-800:]}")
            raise SystemExit(1)
        r_veg, p_veg = float(mm.group(1)), float(mm.group(2))
        prim = pd.read_csv(TAB / f"concordance_mantel_{panel.lower()}.csv")
        r_own = float(prim.loc[prim["comparison"] == "genomic~phenomic", "r"].iloc[0])
        rows.append({"panel": panel, "r_bespoke": r_own, "r_vegan": r_veg,
                     "abs_diff": abs(r_own - r_veg), "p_vegan": p_veg})
        print(f"[33a] {panel}: bespoke r = {r_own:.4f} vs vegan r = {r_veg:.4f} "
              f"(|diff| = {abs(r_own - r_veg):.2e}), vegan P = {p_veg}")

        # ---- (b) Procrustes reconciliation with stage-09 inputs ----
        pca_g = pd.read_csv(TAB / f"pca_{panel.lower()}.csv")
        pca_g["sample"] = pca_g["sample"].astype(str).str.strip()
        G = pca_g.set_index("sample").reindex(common)[
            [c for c in pca_g.columns if c.startswith("PC")]].to_numpy()
        max_pc = G.shape[1]
        P = PCA(n_components=max_pc, random_state=42).fit_transform(Xs)
        keep = ~np.isnan(G).any(axis=1)
        G, P = G[keep], P[keep]
        n = G.shape[0]
        prows = []
        for npc in range(2, max_pc + 1):
            _, _, disp = scipy_procrustes(G[:, :npc], P[:, :npc])
            count = 0
            for _ in range(199):
                perm = RNG.permutation(n)
                _, _, dp_ = scipy_procrustes(G[:, :npc], P[perm, :npc])
                count += dp_ <= disp
            prows.append({"n_pc": npc, "m2": float(disp), "p_perm": (count + 1) / 200,
                          "inputs": "stage09_pipeline (SNPRelate PCs x sklearn feature PCA)"})
        pd.DataFrame(prows).to_csv(TAB / f"robustness_procrustes_{panel.lower()}.csv",
                                   index=False)
        band = [f"{r_['m2']:.3f}" for r_ in prows]
        print(f"[33b] {panel}: M2 (stage-09 inputs) across n_pc 2-{max_pc} = "
              f"{', '.join(band)}")

    pd.DataFrame(rows).to_csv(TAB / "mantel_vegan_validation.csv", index=False)
    print("[33] done")


if __name__ == "__main__":
    main()
