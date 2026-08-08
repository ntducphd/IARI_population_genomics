#!/usr/bin/env python
"""build_source_data.py — Source_Data.xlsx: one worksheet per data-bearing main figure
(house convention). The End matter promises a Source Data
file; this script makes that promise true and keeps it true on rebuild.

Fig 1 (design schematic) has no data sheet. Every other main figure maps to the result tables its
panels are drawn from (per MANIFEST.md); sheets are named FigN[_panel] and carry a first-row
provenance note (source CSV path relative to the compendium root).

Output: analysis/results/source_data/Source_Data.xlsx (copied into the bundle by
the packaging step).
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import ROOT, TAB

OUT_DIR = ROOT / "analysis/results/source_data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT = OUT_DIR / "Source_Data.xlsx"

# sheet name -> source table filename (in analysis/results/tables/)
SHEETS = {
    "Fig2a_pca_set1": "pca_set1.csv",
    "Fig2b_pca_set2": "pca_set2.csv",
    "Fig2cd_admix_cv": ["admixture_set1_cv.csv", "admixture_set2_cv.csv"],
    "Fig2e_Q_set1": "admixture_set1_Q.csv",
    "Fig2f_Q_set2": "admixture_set2_Q.csv",
    "Fig3ab_nj_tips": ["nj_set1_tips.csv", "nj_set2_tips.csv"],
    "Fig4ab_gst_pairwise": ["fst_pairwise_set1.csv", "fst_pairwise_set2.csv"],
    "Fig4c_diversity": "diversity_summary.csv",
    "Fig4d_ld_decay": ["ld_decay_set1.csv", "ld_decay_set2.csv"],
    "Fig5a_mantel_ci": ["concordance_mantel_set1.csv", "concordance_mantel_set2.csv",
                        "mantel_bootstrap_set1.csv", "mantel_bootstrap_set2.csv"],
    "Fig5b_procrustes_band": ["robustness_procrustes_set1.csv",
                              "robustness_procrustes_set2.csv"],
    "Fig5c_classifier_null": ["classifier_baselines_set1.csv", "classifier_baselines_set2.csv",
                              "external_label_classifier_set1.csv"],
    "Fig5d_feature_families": ["concordance_feature_attribution_set1.csv",
                               "concordance_feature_attribution_set2.csv"],
    "Fig5e_stature_collapse": ["confound_mantel_set1.csv", "confound_mantel_set2.csv"],
    "Fig5f_subset_classifiers": ["classifier_baselines_set1.csv",
                                 "classifier_baselines_set2.csv"],
    "Fig6a_roh_classes": ["roh_length_classes_set1.csv", "roh_length_classes_set2.csv"],
    "Fig6b_froh_clusters": "roh_indiv_set1.csv",
    "Fig6c_ne_trajectory": ["ne_trajectory_set1.csv", "ne_trajectory_set2.csv",
                            "ne_interval_summary.csv"],
    "Fig6d_sfs": "sfs_folded_set1.csv",
    "Fig7ab_pcadapt": "pcadapt_summary.csv",   # full per-SNP tables exceed xlsx practicality;
                                                # summary here, per-SNP CSVs in the data deposit
    "Fig7c_ihs": "ihs_set1.csv",
    "Fig7d_xpehh": "xpehh_set1.csv",
    "Fig7_overlap_null": "overlap_null_summary.csv",
    "Fig8a_core_curves": ["core_collection_curve_set1.csv", "core_collection_curve_set2.csv"],
    "Fig8bc_private_richness": ["richness_set1.csv", "richness_set2.csv"],
    "Fig9ab_height_by_cluster": ["height_by_cluster_set1.csv", "height_by_cluster_set2.csv",
                                 "height_by_cluster_summary.csv"],
    "Fig9cd_distance_pairs": ["distance_pairs_set1.csv", "distance_pairs_set2.csv"],
    "Fig9e_external_confusion": ["external_confusion_set1.csv",
                                 "external_label_classifier_set1.csv"],
    "Fig9f_robustness_ladder": ["classifier_baselines_set1.csv",
                                "classifier_baselines_set2.csv",
                                "external_label_classifier_set1.csv"],
    "Fig10a_scoreboard": "replication_scoreboard.csv",
    "Fig10b_maf_spectrum": ["maf_spectrum_set1.csv", "maf_spectrum_set2.csv"],
    "Fig10c_snp_density": "snp_density_by_chrom.csv",
}
MAX_ROWS = 100_000   # xlsx practicality guard; larger tables are truncated with a note


def main():
    with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
        for sheet, src in SHEETS.items():
            srcs = src if isinstance(src, list) else [src]
            frames = []
            for s in srcs:
                p = TAB / s
                if not p.exists():
                    print(f"  (skip {sheet}: {s} missing)")
                    frames = []
                    break
                df = pd.read_csv(p)
                if len(srcs) > 1:
                    df.insert(0, "source_file", s)
                frames.append(df)
            if not frames:
                continue
            df = pd.concat(frames, ignore_index=True)
            note = f"Source: analysis/results/tables/{'; '.join(srcs)}"
            truncated = len(df) > MAX_ROWS
            if truncated:
                df = df.head(MAX_ROWS)
                note += f" (first {MAX_ROWS} rows; full table in the data deposit)"
            meta = pd.DataFrame({df.columns[0]: [note]})
            meta.to_excel(xw, sheet_name=sheet[:31], index=False, header=False)
            df.to_excel(xw, sheet_name=sheet[:31], index=False, startrow=2)
    print(f"[source-data] -> {OUT}")


if __name__ == "__main__":
    main()
