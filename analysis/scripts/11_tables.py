#!/usr/bin/env python
"""11_tables.py — assemble manuscript-ready tables (markdown + csv) from the individual result
CSVs written by Stages 2-9. Main tables go to results/tables/Table_*, supplementary to
results/supp_tables/SuppTable_*.
"""
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import TAB, STAB

PANELS = ["Set1", "Set2"]


def md_table(df, caption):
    return f"**{caption}**\n\n" + df.to_markdown(index=False, floatfmt=".4g") + "\n"


# ---- Table 1: panel summary + diversity ----
prep = pd.read_csv(TAB / "prep_summary.csv")
div = pd.read_csv(TAB / "diversity_summary.csv")
bestk = pd.read_csv(TAB / "admixture_bestK.csv")
t1 = prep.merge(div, on="panel").merge(bestk, on=["panel", "n"], how="left")
t1 = t1[["panel", "platform", "n", "snps_qc", "snps_pruned", "bestK_elbow", "bestK_global_min",
         "mean_PIC", "mean_Ho", "mean_He", "mean_F", "genome_pi", "genome_theta_w", "genome_tajima_d"]]
t1.columns = ["Panel", "Platform", "n", "SNPs (QC)", "SNPs (pruned)", "K (elbow)", "K (global-min)",
              "Mean PIC", "Mean Ho", "Mean He", "Mean F", "Genome pi", "Genome theta_W", "Tajima's D"]
t1.to_csv(TAB / "Table_1_panel_diversity_summary.csv", index=False)
with open(TAB / "Table_1_panel_diversity_summary.md", "w", encoding="utf-8") as f:
    f.write(md_table(t1, "Table 1. Panel summary, admixture K-selection, and genetic diversity."))
print("Table 1 written")

# ---- Table 2: Fst + AMOVA ----
fst_g = pd.concat([pd.read_csv(TAB / f"fst_global_{p.lower()}.csv") for p in PANELS])
amova = pd.concat([pd.read_csv(TAB / f"amova_{p.lower()}.csv") for p in PANELS])
amova_wide = amova.pivot_table(index="panel", columns="stratum", values="pct_variance").reset_index()
t2 = fst_g.merge(amova_wide, on="panel")
t2 = t2[["panel", "n", "n_clusters", "global_fst", "perm_p", "n_perm", "cluster", "Error"]]
t2.columns = ["Panel", "n", "Admixture clusters (K)", "Global Fst", "Fst permutation P", "N permutations",
              "AMOVA % variance among clusters", "AMOVA % variance within clusters"]
t2.to_csv(TAB / "Table_2_fst_amova.csv", index=False)
with open(TAB / "Table_2_fst_amova.md", "w", encoding="utf-8") as f:
    f.write(md_table(t2, "Table 2. Genetic differentiation (Fst, between admixture clusters) and AMOVA variance partition."))
print("Table 2 written")

# ---- Table 3: LD decay + panel platform SNP density ----
ld = pd.read_csv(TAB / "ld_decay_summary.csv")
ld = ld[["panel", "n_snps_used", "thinned", "r2_at_shortest_bin", "half_decay_bp"]].copy()
ld["thinned"] = ld["thinned"].map({True: "Yes (~10%)", False: "No (full QC set)"})
ld["half_decay_bp"] = (ld["half_decay_bp"] / 1000).round(0).astype(int)
ld.columns = ["Panel", "SNPs used", "Thinned before r2", "r2 at shortest bin (12.5 kb)", "Half-decay distance (kb)"]
ld.to_csv(TAB / "Table_3_ld_decay.csv", index=False)
with open(TAB / "Table_3_ld_decay.md", "w", encoding="utf-8") as f:
    f.write(md_table(ld, "Table 3. Linkage-disequilibrium decay."))
print("Table 3 written")

# ---- Table 4: Pillar B concordance (the key results table) ----
mant = pd.concat([pd.read_csv(TAB / f"concordance_mantel_{p.lower()}.csv") for p in PANELS])
# accession-bootstrap 95% CI (stage 30) for the primary genomic~phenomic row of each panel
ci_map = {}
for p in PANELS:
    bp = TAB / f"mantel_bootstrap_{p.lower()}.csv"
    if bp.exists():
        b = pd.read_csv(bp).iloc[0]
        ci_map[p] = f"[{b.ci95_lo:.3f}, {b.ci95_hi:.3f}]"
mant["bootstrap_95CI"] = [
    ci_map.get(r.panel, "--") if r.comparison == "genomic~phenomic" else "--"
    for r in mant.itertuples()]
mant.to_csv(TAB / "Table_4a_mantel_concordance.csv", index=False)
mant_disp = mant.copy()
mant_disp.columns = ["Panel", "n", "Comparison", "Mantel r", "P", "Bootstrap 95% CI"]
attr = pd.concat([pd.read_csv(TAB / f"concordance_feature_attribution_{p.lower()}.csv") for p in PANELS])
attr.to_csv(TAB / "Table_4b_feature_attribution.csv", index=False)
attr_disp = attr.copy()
attr_disp.columns = ["Panel", "Feature family", "N features", "Genomic Mantel r", "P"]
cls = pd.concat([pd.read_csv(TAB / f"concordance_classification_{p.lower()}.csv") for p in PANELS
                  if (TAB / f"concordance_classification_{p.lower()}.csv").exists()])
cls.to_csv(TAB / "Table_4c_classification.csv", index=False)
cls_disp = cls.copy()
cls_disp.columns = ["Panel", "n", "Admixture clusters (K)", "CV folds", "RF CV accuracy", "Majority-class baseline"]
proc = pd.concat([pd.read_csv(TAB / f"concordance_procrustes_{p.lower()}.csv") for p in PANELS])
proc.to_csv(TAB / "Table_4d_procrustes.csv", index=False)
proc_disp = proc.copy()
proc_disp.columns = ["Panel", "n", "Procrustes M2", "Permutation P", "N permutations"]
conf = pd.concat([pd.read_csv(TAB / f"concordance_confounder_{p.lower()}.csv") for p in PANELS])
conf.to_csv(TAB / "Table_4e_confounder.csv", index=False)
conf_disp = conf.copy()
conf_disp.columns = ["Panel", "Raw phenomic~trait r", "Raw P", "Partial r (| genomic)", "Partial P", "% change"]
with open(TAB / "Table_4_pillarB_concordance.md", "w", encoding="utf-8") as f:
    f.write(md_table(mant_disp, "Table 4a. Mantel and partial-Mantel tests: genomic, phenomic, and trait distances."))
    f.write("\n" + md_table(proc_disp, "Table 4b. Procrustes/PROTEST: genomic-PCA vs phenomic-PCA."))
    f.write("\n" + md_table(cls_disp, "Table 4c. Supervised classification: phenomic features -> admixture cluster."))
    f.write("\n" + md_table(attr_disp, "Table 4d. Feature-type attribution: which imaging channel carries genomic signal."))
    f.write("\n" + md_table(conf_disp, "Table 4e. Structure-as-confounder test: phenomic~trait, raw vs partialled for genomic distance."))
print("Table 4 (a-e) written")

# ---- Table 5: core collection summary ----
core = pd.concat([pd.read_csv(TAB / f"core_collection_summary_{p.lower()}.csv") for p in PANELS])
core = core[["panel", "n_total", "recommended_core_size", "pct_of_panel", "pct_diversity_retained"]]
core.columns = ["Panel", "n (total)", "Recommended core size", "% of panel", "% diversity retained"]
core.to_csv(TAB / "Table_5_core_collection.csv", index=False)
with open(TAB / "Table_5_core_collection.md", "w", encoding="utf-8") as f:
    f.write(md_table(core, "Table 5. Recommended core collection (M-strategy, max-min diversity)."))
print("Table 5 written")

# ==== Supplementary tables ====
# S1: full panel prep detail
prep.to_csv(STAB / "SuppTable_S1_prep_summary.csv", index=False)

# S2: subpopulation capture (Set1, ground-truth 3K-RGP labels)
cap = pd.read_csv(TAB / "subpop_capture_set1.csv")
cap.to_csv(STAB / "SuppTable_S2_subpop_capture_set1.csv", index=False)

# S3: full pairwise Fst
pw = pd.concat([pd.read_csv(TAB / f"fst_pairwise_{p.lower()}.csv") for p in PANELS])
pw.to_csv(STAB / "SuppTable_S3_pairwise_fst.csv", index=False)

# S4: MAF spectrum
maf = pd.concat([pd.read_csv(TAB / f"maf_spectrum_{p.lower()}.csv") for p in PANELS])
maf.to_csv(STAB / "SuppTable_S4_maf_spectrum.csv", index=False)

# S5: structure consensus (ARI/NMI, cross-method)
sc = pd.read_csv(TAB / "structure_consensus_summary.csv")
sc.to_csv(STAB / "SuppTable_S5_structure_consensus.csv", index=False)

# S6: core collection curve (full retention curve, both panels)
curve = pd.concat([pd.read_csv(TAB / f"core_collection_curve_{p.lower()}.csv") for p in PANELS])
curve.to_csv(STAB / "SuppTable_S6_core_collection_curve.csv", index=False)

# S7: underused accessions (pre-breeding priority list)
under = pd.concat([pd.read_csv(TAB / f"underused_accessions_{p.lower()}.csv") for p in PANELS])
under.to_csv(STAB / "SuppTable_S7_underused_accessions.csv", index=False)

print("Supplementary tables S1-S7 written")

# ============================================================================================
# Flagship additions (stages 13-26, 2026-08-08) — upgrades to Tables 1-3, new Tables 6-7,
# supplementary S8-S12. Each block skips gracefully (with a WARNING) if its stage has not run.
# ============================================================================================

def _have(*names):
    missing = [n for n in names if not (TAB / n).exists()]
    if missing:
        print(f"  WARNING: skipped (missing {missing})")
        return False
    return True

# ---- Table 1 UPGRADE: no-MAF diversity + F_ROH + selfing rate ----
if _have("diversity_nomaf_summary.csv", "roh_summary.csv", "hwe_summary.csv"):
    nomaf = pd.read_csv(TAB / "diversity_nomaf_summary.csv")
    roh = pd.read_csv(TAB / "roh_summary.csv")[["panel", "mean_f_roh"]]
    hwe = pd.read_csv(TAB / "hwe_summary.csv")[["panel", "selfing_rate_from_mean_F"]]
    t1v2 = (t1.rename(columns={"Panel": "panel"})
              .drop(columns=["Genome pi", "Genome theta_W", "Tajima's D"])
              .merge(nomaf[["panel", "n_snps_nomaf", "genome_pi", "genome_theta_w",
                            "genome_tajima_d"]], on="panel")
              .merge(roh, on="panel").merge(hwe, on="panel"))
    t1v2.columns = list(t1v2.columns[:-6]) + [
        "SNPs (missingness-only)", "Genome pi (no-MAF)", "Genome theta_W (no-MAF)",
        "Tajima's D (no-MAF)", "Mean F_ROH", "Selfing rate s"]
    t1v2 = t1v2.rename(columns={"panel": "Panel"})
    t1v2.to_csv(TAB / "Table_1_panel_diversity_summary.csv", index=False)
    with open(TAB / "Table_1_panel_diversity_summary.md", "w", encoding="utf-8") as f:
        f.write(md_table(t1v2, "Table 1. Panel summary, admixture K-selection, genetic "
                         "diversity (missingness-only-filtered SNPs; SFS-shape statistics for "
                         "Set 2 are ascertainment-biased and platform-relative), genomic "
                         "inbreeding (F_ROH), and estimated selfing rate."))
    print("Table 1 upgraded (no-MAF diversity + F_ROH + selfing)")

# ---- Table 2 UPGRADE: Weir-Cockerham theta leads; Gst secondary; AMOVA + two-level AMOVA ----
if _have("fst_wc_global_set1.csv", "fst_wc_global_set2.csv"):
    wc = pd.concat([pd.read_csv(TAB / f"fst_wc_global_{p.lower()}.csv").assign(panel=p)
                    for p in PANELS])
    wc["WC theta (95% CI)"] = wc.apply(
        lambda r: f"{r['theta_global']:.3f} ({r['ci95_lo']:.3f}-{r['ci95_hi']:.3f})", axis=1)
    t2v2 = t2.rename(columns={"Panel": "panel"}).merge(
        wc[["panel", "WC theta (95% CI)"]], on="panel").rename(columns={"panel": "Panel"})
    t2v2 = t2v2.rename(columns={"Global Fst": "Nei Gst (secondary)"})
    cols = ["Panel", "n", "Admixture clusters (K)", "WC theta (95% CI)",
            "Nei Gst (secondary)", "Fst permutation P",
            "AMOVA % variance among clusters", "AMOVA % variance within clusters"]
    t2v2 = t2v2[[c for c in cols if c in t2v2.columns]]
    t2v2.to_csv(TAB / "Table_2_fst_amova.csv", index=False)
    with open(TAB / "Table_2_fst_amova.md", "w", encoding="utf-8") as f:
        f.write(md_table(t2v2, "Table 2. Genetic differentiation among admixture clusters: "
                         "Weir-Cockerham theta (block-bootstrap 95% CI over 1-Mb blocks; "
                         "selfing inflates all Fst-family estimators, see Methods), Nei Gst, "
                         "and AMOVA variance partition."))
        if (TAB / "amova2_set1.csv").exists():
            a2 = pd.read_csv(TAB / "amova2_set1.csv")
            a2.columns = ["Component", "sigma^2", "% variance", "P"]
            a2["Component"] = ["Among macro-groups (Indica/Japonica/Aus/Admixed)",
                               "Among 3K-RGP subpopulations within macro-group",
                               "Within subpopulations"][:len(a2)]
            # the within-subpopulation (Error) component has no permutation test -> em dash,
            # not a bare "nan" in the rendered table; pegas returns 0.000 when the
            # observed value beat all 999 permutations -> report the resolution bound, never
            # a literal impossible P = 0
            a2["P"] = a2["P"].map(lambda v: "-" if pd.isna(v)
                                  else ("<0.001" if v < 0.001 else f"{v:.3f}"))
            f.write("\n" + md_table(a2, "Table 2b. Two-level hierarchical AMOVA (Set 1, "
                                    "external 3K-RGP subpopulation labels)."))
    print("Table 2 upgraded (WC theta + two-level AMOVA)")

# ---- Table 3 UPGRADE: standard-threshold decay + empirical thinning check ----
if _have("ld_threshold_summary.csv", "ld_thinning_check.csv"):
    thr = pd.read_csv(TAB / "ld_threshold_summary.csv")
    ld2 = pd.read_csv(TAB / "ld_decay_summary.csv")
    t3v2 = ld2.merge(thr[["panel", "dist_r2_0.2_bp", "dist_r2_0.1_bp"]], on="panel")
    t3v2["half_decay_kb"] = (t3v2["half_decay_bp"] / 1000).round(0).astype(int)
    # a curve that never crosses the threshold within the 1-Mb window is "not reached",
    # never a rendered "nan" (Set 2 mean r2 stays >0.1 through 1 Mb)
    for c in ["dist_r2_0.2_bp", "dist_r2_0.1_bp"]:
        t3v2[c] = t3v2[c].map(lambda v: "not reached (>1,000)" if pd.isna(v)
                              else f"{v / 1000:.0f}")
    t3v2["thinned"] = t3v2["thinned"].map({True: "Yes (~10%)", False: "No (full QC set)"})
    t3v2 = t3v2[["panel", "n_snps_used", "thinned", "r2_at_shortest_bin", "half_decay_kb",
                 "dist_r2_0.2_bp", "dist_r2_0.1_bp"]]
    t3v2.columns = ["Panel", "SNPs used", "Thinned before r2", "r2 at shortest bin",
                    "Half-decay (kb)", "Distance to r2=0.2 (kb)", "Distance to r2=0.1 (kb)"]
    t3v2.to_csv(TAB / "Table_3_ld_decay.csv", index=False)
    chk = pd.read_csv(TAB / "ld_thinning_check.csv")
    chk_disp = chk[["variant", "half_decay_bp", "dist_r2_0.2_bp"]].copy()
    for c in ["half_decay_bp", "dist_r2_0.2_bp"]:
        chk_disp[c] = (chk_disp[c] / 1000).round(0)
    chk_disp.columns = ["Variant (chr 1, 1-10 Mb)", "Half-decay (kb)", "Distance to r2=0.2 (kb)"]
    with open(TAB / "Table_3_ld_decay.md", "w", encoding="utf-8") as f:
        f.write(md_table(t3v2, "Table 3. Linkage-disequilibrium decay, with standard-threshold "
                         "crossing distances."))
        f.write("\n" + md_table(chk_disp, "Table 3b. Empirical thinning-bias check: 10% marker "
                                "thinning over-estimates the Set 1 half-decay distance by ~14%."))
    print("Table 3 upgraded (thresholds + thinning check)")

# ---- Table 6 (NEW): genome history — ROH + Ne INTERVAL ----
if _have("roh_summary.csv", "ne_interval_summary.csv"):
    roh_full = pd.read_csv(TAB / "roh_summary.csv")
    ni = pd.read_csv(TAB / "ne_interval_summary.csv")
    agg = ni.groupby("panel").agg(
        pan_lo=("recent_ne_panmictic", "min"), pan_hi=("recent_ne_panmictic", "max"),
        adj_lo=("recent_ne_selfing_adj", "min"), adj_hi=("recent_ne_selfing_adj", "max"),
        epoch=("epoch_generations_selfing_adj", "mean")).reset_index()
    agg["Recent Ne (panmictic)"] = agg.apply(
        lambda r: f"{r['pan_lo']:.0f}-{r['pan_hi']:.0f}", axis=1)
    agg["Recent Ne (selfing-adj.)"] = agg.apply(
        lambda r: f"{r['adj_lo']:.0f}-{r['adj_hi']:.0f}", axis=1)
    agg["Epoch (generations, selfing-adj.)"] = agg["epoch"].round(0).astype(int)
    t6 = roh_full[["panel", "n_segments", "mean_f_roh", "median_f_roh",
                   "mean_total_roh_mb"]].merge(
        agg[["panel", "Recent Ne (panmictic)", "Recent Ne (selfing-adj.)",
             "Epoch (generations, selfing-adj.)"]], on="panel")
    t6.columns = ["Panel", "ROH segments", "Mean F_ROH", "Median F_ROH",
                  "Mean total ROH (Mb)", "Recent Ne (panmictic)",
                  "Recent Ne (selfing-adj.)", "Epoch (generations, selfing-adj.)"]
    t6.to_csv(TAB / "Table_6_genome_history.csv", index=False)
    with open(TAB / "Table_6_genome_history.md", "w", encoding="utf-8") as f:
        f.write(md_table(t6, "Table 6. Genome history: runs of homozygosity and "
                         "linkage-disequilibrium-based recent effective population size, "
                         "reported as intervals across the inbreeding-input (PLINK F vs "
                         "F_ROH) and map-density (3-5 cM/Mb) grid; full grid in Supp. Table "
                         "S15; assumption ledger in Methods."))
    print("Table 6 written (genome history, Ne intervals)")

# ---- Table 7 (NEW): selection-scan summary ----
if _have("pcadapt_summary.csv"):
    pc7 = pd.read_csv(TAB / "pcadapt_summary.csv")
    pc7 = pc7.rename(columns={"panel": "Panel", "K_used": "pcadapt K", "n_snps": "SNPs",
                              "n_outliers": "pcadapt outliers (BH q<0.05)",
                              "prop_outliers": "Proportion"})
    with open(TAB / "Table_7_selection.md", "w", encoding="utf-8") as f:
        f.write(md_table(pc7, "Table 7a. Fst-outlier selection scan (pcadapt)."))
        if (TAB / "haplo_scan_summary.csv").exists():
            hs = pd.read_csv(TAB / "haplo_scan_summary.csv")
            # reader-facing headers (raw code-style names leaked into the built DOCX -- QC catch)
            hs = hs.rename(columns={
                "panel": "Panel", "n_snps_ihs": "iHS SNPs", "n_outlier_ihs": "iHS outliers",
                "n_snps_xpehh": "XP-EHH SNPs", "n_outlier_xpehh": "XP-EHH outliers",
                "popA": "Group A", "nA": "n (A)", "popB": "Group B", "nB": "n (B)",
                "n_regions": "Candidate regions"})
            f.write("\n" + md_table(hs, "Table 7b. Haplotype-based selection scans (Set 1, "
                                    "Beagle-phased; iHS panel-wide, XP-EHH between the two "
                                    "largest admixture clusters; |score| >= 3 outliers; "
                                    "selfing caveat in Methods)."))
        else:
            print("  (Table 7b pending stage 17)")
        if (TAB / "overlap_null_summary.csv").exists():
            on = pd.read_csv(TAB / "overlap_null_summary.csv")
            on7 = on[["observed_fraction", "null_mean_fraction", "empirical_p", "enrichment",
                      "exact_overlap_fraction", "pcadapt_outlier_density",
                      "exact_snp_enrichment"]].copy()
            on7.columns = ["Observed region overlap", "Rotation-null mean", "Empirical P",
                           "Enrichment", "Exact-SNP overlap", "pcadapt outlier density",
                           "Exact-SNP enrichment"]
            f.write("\n" + md_table(on7, "Table 7c. Chance expectation for two-scan "
                                    "convergence (999 within-chromosome rotations): the "
                                    "region-level overlap is chance-level and the exact-SNP "
                                    "overlap is depleted — the two scans are complementary "
                                    "screens, not cross-validated hits."))
    pc7.to_csv(TAB / "Table_7_selection.csv", index=False)
    print("Table 7 written (selection)")

# ---- Supplementary S8-S12 ----
if _have("f3_set1.csv"):
    pd.read_csv(TAB / "f3_set1.csv").to_csv(STAB / "SuppTable_S8_f3_admixture_tests.csv",
                                            index=False)
if _have("richness_set1.csv", "richness_set2.csv"):
    pd.concat([pd.read_csv(TAB / f"richness_{p.lower()}.csv").assign(panel=p)
               for p in PANELS]).to_csv(STAB / "SuppTable_S9_allelic_richness.csv", index=False)
if _have("roh_cluster_summary_set1.csv", "roh_cluster_summary_set2.csv"):
    pd.concat([pd.read_csv(TAB / f"roh_cluster_summary_{p.lower()}.csv").assign(panel=p)
               for p in PANELS]).to_csv(STAB / "SuppTable_S10_roh_by_cluster.csv", index=False)
if _have("robustness_summary.csv", "robustness_fdr_mantel.csv"):
    pd.read_csv(TAB / "robustness_summary.csv").to_csv(
        STAB / "SuppTable_S11_robustness_summary.csv", index=False)
    pd.read_csv(TAB / "robustness_fdr_mantel.csv").to_csv(
        STAB / "SuppTable_S12_fdr_ledger.csv", index=False)

# ---- S13-S16: review-response stages 30-33 (2026-08-08) ----
if _have("confound_mantel_set1.csv", "confound_mantel_set2.csv"):
    pd.concat([pd.read_csv(TAB / f"confound_mantel_{p.lower()}.csv") for p in PANELS] +
              [pd.read_csv(TAB / f"spatial_proxy_{p.lower()}.csv") for p in PANELS
               if (TAB / f"spatial_proxy_{p.lower()}.csv").exists()]).to_csv(
        STAB / "SuppTable_S13_confound_tests.csv", index=False)
if _have("classifier_baselines_set1.csv", "classifier_baselines_set2.csv"):
    frames = [pd.read_csv(TAB / f"classifier_baselines_{p.lower()}.csv") for p in PANELS]
    if (TAB / "external_label_classifier_set1.csv").exists():
        frames.append(pd.read_csv(TAB / "external_label_classifier_set1.csv"))
    pd.concat(frames).to_csv(STAB / "SuppTable_S14_classifier_robustness.csv", index=False)
if _have("ne_interval_summary.csv", "within_cluster_history_set1.csv"):
    pd.read_csv(TAB / "ne_interval_summary.csv").to_csv(
        STAB / "SuppTable_S15_ne_interval.csv", index=False)
if _have("overlap_null_summary.csv", "mantel_vegan_validation.csv"):
    pd.read_csv(TAB / "overlap_null_summary.csv").to_csv(
        STAB / "SuppTable_S16_overlap_null_validation.csv", index=False)
# ---- mechanism/robustness supp tables (stages 37/38) ----
if _have("pst_fst_sweep.csv"):
    pd.read_csv(TAB / "pst_fst_sweep.csv").to_csv(
        STAB / "SuppTable_S17_pst_fst_sweep.csv", index=False)
if _have("temporal_mantel.csv"):
    pd.read_csv(TAB / "temporal_mantel.csv").to_csv(
        STAB / "SuppTable_S18_temporal_mantel.csv", index=False)
if _have("classifier_uncertainty.csv"):
    pd.read_csv(TAB / "classifier_uncertainty.csv").to_csv(
        STAB / "SuppTable_S19_classifier_uncertainty.csv", index=False)
if _have("sizeleakage_summary.csv", "residualised_classifier.csv", "kbattery_forensics.csv"):
    blocks = []
    for src in ["sizeleakage_summary.csv", "residualised_classifier.csv",
                "kbattery_forensics.csv"]:
        b = pd.read_csv(TAB / src)
        b.insert(0, "block", src.replace(".csv", ""))
        blocks.append(b)
    pd.concat(blocks, ignore_index=True).to_csv(
        STAB / "SuppTable_S20_mechanism_partition.csv", index=False)

print(f"\n-> {TAB} (Table 1-7) + {STAB} (SuppTable S1-S20)")
