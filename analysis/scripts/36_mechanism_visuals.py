#!/usr/bin/env python
"""36_mechanism_visuals.py — [stage 36] source tables for the two NEW main figures added in the
final 10-figure scheme (author direction 2026-08-08: genuinely new figures, not supp promotions).

Fig 9 "Mechanism made visible" needs:
  height_by_cluster_{set}.csv   per-accession height index (mean z-score of the 9 height-derived
                                features, stage-30 definition) + admixture cluster; Kruskal-Wallis
                                P across clusters in height_by_cluster_summary.csv
  distance_pairs_{set}.csv      upper-triangle (genomic IBS distance, phenomic Euclidean distance)
                                pairs — the raw material behind the primary Mantel r
  external_confusion_set1.csv   confusion matrix of the stage-30 external-label classifier
                                (identical protocol + seeds; accuracy is asserted to match the
                                stage-30 headline 58.9% so the figure can never drift from Table S14)

Fig 10 "Two platforms, one biology" needs:
  replication_scoreboard.csv    the paper's dimensionless headline statistics side by side
                                (Set 1 vs Set 2), assembled from the existing summary tables —
                                nothing recomputed, chỉ gom số đã có (không bịa)

Reuses stage-30 helpers via importlib (module name starts with a digit).
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal
from sklearn.metrics import confusion_matrix

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths import TAB

_spec = importlib.util.spec_from_file_location("stage30", HERE / "30_confound_robustness.py")
s30 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s30)

import re


def height_cols_of(coh):
    img_cols = [c for c in coh.columns if c.startswith("img")]
    return [c for c in img_cols if ("Height" in c) or re.search(r"_iPH\b", c)], img_cols


def main():
    kw_rows = []
    for panel in ["Set1", "Set2"]:
        common, coh, Dg, y, _ = s30.load_panel(panel)
        height_cols, img_cols = height_cols_of(coh)

        # ---- Fig 9a,b: height index (mean z-score over the stage-30 height features) ----
        H = coh[height_cols].to_numpy(dtype=float)
        H = np.nan_to_num(H, nan=np.nanmean(H))
        Hz = (H - H.mean(axis=0)) / H.std(axis=0)
        hidx = Hz.mean(axis=1)
        pd.DataFrame({"sample": common, "cluster": y,
                      "height_index_z": hidx}).to_csv(
            TAB / f"height_by_cluster_{panel.lower()}.csv", index=False)
        groups = [hidx[y == c] for c in np.unique(y) if (y == c).sum() >= 2]
        kw_stat, kw_p = kruskal(*groups)
        kw_rows.append({"panel": panel, "n": len(common),
                        "n_height_features": len(height_cols),
                        "kruskal_H": float(kw_stat), "kruskal_p": float(kw_p)})
        print(f"[36] {panel}: height index by cluster, Kruskal-Wallis H = {kw_stat:.1f}, "
              f"P = {kw_p:.2e}")

        # ---- Fig 9c,d: the raw distance-pair cloud behind the primary Mantel r ----
        Dp = s30.dist_std(coh[img_cols].to_numpy(dtype=float))
        iu = np.triu_indices_from(Dg, k=1)
        pd.DataFrame({"genomic_dist": Dg[iu], "phenomic_dist": Dp[iu]}).to_csv(
            TAB / f"distance_pairs_{panel.lower()}.csv", index=False)
        print(f"[36] {panel}: {len(iu[0]):,} distance pairs written")

        # ---- Fig 9e: external-label confusion matrix (Set 1; stage-30 protocol verbatim) ----
        if panel == "Set1":
            from sklearn.preprocessing import StandardScaler
            X_all = coh[img_cols].to_numpy(dtype=float)
            Xs_all = StandardScaler().fit_transform(
                np.nan_to_num(X_all, nan=np.nanmean(X_all)))
            sub = pd.read_csv(TAB / "subpop_assignment_set1.csv")
            m = dict(zip(sub["sample"].astype(str), sub["subpopulation"].astype(str)))
            lab = pd.Series([m.get(s, "NA") for s in common])
            counts = lab.value_counts()
            keepg = lab.isin([g for g in counts.index if g != "NA" and counts[g] >= 5])
            ye, names = pd.factorize(lab[keepg])
            Xe = Xs_all[keepg.to_numpy()]
            acc_e, ke, _, _, pred = s30.clf_acc(Xe, ye)
            headline = pd.read_csv(TAB / "external_label_classifier_set1.csv").iloc[0]
            assert abs(acc_e - float(headline.accuracy)) < 1e-9, (
                f"external confusion acc {acc_e} != stage-30 headline {headline.accuracy}")
            cm = confusion_matrix(ye, pred)
            rows = [{"true_group": names[i], "pred_group": names[j], "count": int(cm[i, j])}
                    for i in range(len(names)) for j in range(len(names))]
            pd.DataFrame(rows).to_csv(TAB / "external_confusion_set1.csv", index=False)
            print(f"[36] Set1 external confusion written (acc = {acc_e:.3f} == stage-30, "
                  f"{ke}-fold CV, groups: {', '.join(names)})")

    pd.DataFrame(kw_rows).to_csv(TAB / "height_by_cluster_summary.csv", index=False)

    # ---- Fig 10a: replication scoreboard (dimensionless headline stats, both panels) ----
    div = pd.read_csv(TAB / "diversity_summary.csv").set_index("panel")
    hwe = pd.read_csv(TAB / "hwe_summary.csv").set_index("panel")
    roh = pd.read_csv(TAB / "roh_summary.csv").set_index("panel")
    rows = []

    def add(metric, s1, s2, family):
        rows.append({"metric": metric, "set1": float(s1), "set2": float(s2),
                     "family": family})

    add("Observed heterozygosity (Ho)", div.loc["Set1", "mean_Ho"],
        div.loc["Set2", "mean_Ho"], "Diversity & mating")
    add("Expected heterozygosity (He)", div.loc["Set1", "mean_He"],
        div.loc["Set2", "mean_He"], "Diversity & mating")
    add("Inbreeding coefficient (F)", div.loc["Set1", "mean_F"],
        div.loc["Set2", "mean_F"], "Diversity & mating")
    add("Selfing rate (s)", hwe.loc["Set1", "selfing_rate_from_mean_F"],
        hwe.loc["Set2", "selfing_rate_from_mean_F"], "Diversity & mating")
    add("Genomic inbreeding (F_ROH)", roh.loc["Set1", "mean_f_roh"],
        roh.loc["Set2", "mean_f_roh"], "Genome history")
    g1 = pd.read_csv(TAB / "fst_global_set1.csv").iloc[0]
    g2 = pd.read_csv(TAB / "fst_global_set2.csv").iloc[0]
    add("Nei Gst (between clusters)", g1.global_fst, g2.global_fst, "Differentiation")
    w1 = pd.read_csv(TAB / "fst_wc_global_set1.csv").iloc[0]
    w2 = pd.read_csv(TAB / "fst_wc_global_set2.csv").iloc[0]
    add("Weir-Cockerham theta", w1.theta_global, w2.theta_global, "Differentiation")
    m1 = pd.read_csv(TAB / "concordance_mantel_set1.csv")
    m2 = pd.read_csv(TAB / "concordance_mantel_set2.csv")
    add("Mantel r (genomic~phenomic)",
        m1.loc[m1.comparison == "genomic~phenomic", "r"].iloc[0],
        m2.loc[m2.comparison == "genomic~phenomic", "r"].iloc[0], "Concordance")
    c1 = pd.read_csv(TAB / "classifier_baselines_set1.csv")
    c2 = pd.read_csv(TAB / "classifier_baselines_set2.csv")
    sel = lambda df: df[(df.features == "full_204") &
                        (df.target == "admixture_elbowK")].iloc[0]
    add("Classifier accuracy (imaging -> cluster)", sel(c1).accuracy, sel(c2).accuracy,
        "Concordance")
    add("Classifier permutation null", sel(c1).null_mean, sel(c2).null_mean, "Concordance")
    a1 = pd.read_csv(TAB / "structure_consensus_set1.csv")
    a2 = pd.read_csv(TAB / "structure_consensus_set2.csv")
    ari = lambda df: df[(df.method_a.str.startswith("UMAP")) &
                        (df.method_b == "admixture_argmaxQ")]["ARI"].iloc[0]
    add("ARI (UMAP vs admixture)", ari(a1), ari(a2), "Structure")
    k1 = pd.read_csv(TAB / "core_collection_summary_set1.csv").iloc[0]
    k2 = pd.read_csv(TAB / "core_collection_summary_set2.csv").iloc[0]
    add("10% core diversity retained (/100)", k1.pct_diversity_retained / 100,
        k2.pct_diversity_retained / 100, "Application")
    pd.DataFrame(rows).to_csv(TAB / "replication_scoreboard.csv", index=False)
    print(f"[36] replication scoreboard: {len(rows)} metrics written")
    print("[36] done")


if __name__ == "__main__":
    main()
