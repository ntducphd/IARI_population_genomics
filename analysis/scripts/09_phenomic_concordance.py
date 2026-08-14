#!/usr/bin/env python
"""09_phenomic_concordance.py — PILLAR B: does canopy-imaging phenomic data recover genetic
population structure? This is the primary, broad-audience result of the paper
(Pillar B is the paper's hook, not Pillar A/C).

Implemented directly in Python rather than via R's vegan: Mantel/partial-Mantel and Procrustes/
PROTEST are simple, well-defined textbook procedures (permutation tests on distance/configuration
matrices), and this session's environment audit found this R 4.4.3 build unreliable even for much
simpler operations (see hierfstat, adegenet findings in the environment-audit record) -- a from-scratch,
easily-verified Python implementation is more trustworthy here than fighting vegan's exact API
under a flaky interpreter for the single most important analysis in the paper.

Three distances/configurations per panel:
  genomic  = 1 - IBS (from Stage 3, 04_kinship_tree.R), restricted to phenotyped accessions
  phenomic = Euclidean distance on standardised 204 canopy-imaging features (cohort_set{1,2}.csv)
  trait    = Euclidean distance on standardised traditional NUE traits (same cohort file)

Analyses:
  1. Mantel: genomic<->phenomic, genomic<->trait (Pearson r on upper-triangle, permutation p)
  2. Partial Mantel: genomic<->phenomic | trait, and the reverse
  3. Procrustes + PROTEST: genomic PCA (Stage 1 PC1-4) vs phenomic PCA (PC1-4 of the 204 features)
  4. Supervised classification: RandomForest, phenomic features -> admixture cluster label
     (Stage 2's elbow-K argmax-Q clusters), 5-fold stratified CV accuracy vs majority-class baseline
  5. Feature-type attribution: repeat the genomic<->phenomic Mantel test using only Colour, only
     NIR, or only Size/morphology features (three feature families identified from the 204 columns)
  6. Structure-as-confounder: Mantel(phenomic, trait) vs partial Mantel(phenomic, trait | genomic)
     -- if the partial correlation drops sharply, genomic structure is inflating the raw
     phenomic-trait relationship that the companion phenomic-selection study's accuracy claims rest on

Outputs (analysis/results/tables/): concordance_mantel_{set}.csv, concordance_procrustes_{set}.csv,
  concordance_classification_{set}.csv, concordance_feature_attribution_{set}.csv,
  concordance_confounder_{set}.csv
"""
import sys
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import accuracy_score
from sklearn.pipeline import make_pipeline

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import TAB, COHORT_SET1, COHORT_SET2

SEED = 42
N_PERM = 999
rng = np.random.default_rng(SEED)

COHORT = {"Set1": COHORT_SET1, "Set2": COHORT_SET2}
NUE_TRAITS = ["NUEb", "NUpE", "NUtE", "PNUE", "NHI", "NUE1963"]


def upper(d):
    return d[np.triu_indices_from(d, k=1)]


def mantel(d1, d2, n_perm=N_PERM):
    u1, u2 = upper(d1), upper(d2)
    r_obs, _ = pearsonr(u1, u2)
    n = d1.shape[0]
    idx = np.arange(n)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(idx)
        d2p = d2[np.ix_(perm, perm)]
        r_p, _ = pearsonr(u1, upper(d2p))
        if abs(r_p) >= abs(r_obs):
            count += 1
    p = (count + 1) / (n_perm + 1)
    return r_obs, p


def partial_mantel(d1, d2, d3, n_perm=N_PERM):
    """Partial Mantel of d1~d2 controlling for d3: partial correlation of the upper-triangles,
    permutation p-value via permuting d2's rows/cols (the standard partial-Mantel permutation
    scheme, e.g. Smouse et al. 1986)."""
    u1, u2, u3 = upper(d1), upper(d2), upper(d3)
    r12, _ = pearsonr(u1, u2); r13, _ = pearsonr(u1, u3); r23, _ = pearsonr(u2, u3)
    def partial_r(r12, r13, r23):
        denom = np.sqrt((1 - r13**2) * (1 - r23**2))
        return (r12 - r13 * r23) / denom if denom > 0 else np.nan
    r_obs = partial_r(r12, r13, r23)
    n = d1.shape[0]
    idx = np.arange(n)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(idx)
        d2p = d2[np.ix_(perm, perm)]
        u2p = upper(d2p)
        r12p, _ = pearsonr(u1, u2p); r23p, _ = pearsonr(u2p, u3)
        r_p = partial_r(r12p, r13, r23p)
        if np.isfinite(r_p) and abs(r_p) >= abs(r_obs):
            count += 1
    p = (count + 1) / (n_perm + 1)
    return r_obs, p


def protest(x1, x2, n_perm=N_PERM):
    m1, m2, disparity_obs = procrustes(x1, x2)
    n = x1.shape[0]
    idx = np.arange(n)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(idx)
        _, _, d_p = procrustes(x1, x2[perm])
        if d_p <= disparity_obs:
            count += 1
    p = (count + 1) / (n_perm + 1)
    return disparity_obs, p


FEATURE_GROUPS = {
    "Colour": lambda c: "Color." in c,
    "NIR": lambda c: c.startswith(("img68_NIR", "img75_NIR", "img83_NIR",
                                    "img68_IR.", "img75_IR.", "img83_IR.")),
}
def feature_group_of(col):
    for name, test in FEATURE_GROUPS.items():
        if test(col):
            return name
    return "Size_morphology"


def normalize_id(x):
    """Genotype<->phenotype bridge rule (data/README.md): Set2 IDs match CASE-INSENSITIVELY --
    numeric IRIS IDs are consistent either way, but ~38 Set2 accessions are UPPERCASE variety names
    in the genotype file vs Title-Case in the phenotype/cohort file. Normalise by upper-casing and
    stripping spaces/dashes on BOTH sides before joining, exactly as documented, instead of a naive
    exact-string match (which recovers only 109/147 for Set2, the known-incomplete number the
    project's own data/README.md explicitly warns about)."""
    return re.sub(r"[\s\-]", "", str(x).upper())


for panel in ["Set1", "Set2"]:
    print(f"\n{'='*20} [{panel}] Pillar B concordance {'='*20}")
    cohort = pd.read_csv(COHORT[panel])
    cohort["Taxa"] = cohort["Taxa"].astype(str).str.strip()
    cohort["_key"] = cohort["Taxa"].map(normalize_id)

    ibs = pd.read_csv(TAB / f"ibs_dist_{panel.lower()}.csv")
    genomic_ids = sorted(set(ibs["sample_a"]) | set(ibs["sample_b"]))
    genomic_key_to_id = {normalize_id(g): g for g in genomic_ids}

    common_keys = sorted(set(cohort["_key"]) & set(genomic_key_to_id))
    common = [genomic_key_to_id[k] for k in common_keys]   # canonical genomic-side IDs, for indexing
    print(f"  cohort n={len(cohort)}, genomic n={len(genomic_ids)}, common={len(common)} "
          f"(case-insensitive bridge match, per data/README.md)")
    if len(common) < 20:
        print(f"  SKIP {panel}: too few common accessions for a meaningful concordance test")
        continue

    key_to_taxa = dict(zip(cohort["_key"], cohort["Taxa"]))
    coh = cohort.drop_duplicates("_key").set_index("_key").loc[common_keys]
    coh.index = common   # re-key to the canonical genomic-side ID so all downstream lookups align
    n = len(common)

    # ---- genomic distance matrix, restricted + reordered to `common` ----
    ibs_p = ibs.pivot(index="sample_a", columns="sample_b", values="distance")
    ibs_p = ibs_p.reindex(index=common, columns=common)
    Dg = ibs_p.to_numpy().copy()
    np.fill_diagonal(Dg, 0.0)

    # ---- phenomic distance (standardised 204 img features) ----
    img_cols = [c for c in coh.columns if c.startswith("img")]
    Ximg = coh[img_cols].to_numpy(dtype=float)
    Ximg = np.nan_to_num(Ximg, nan=np.nanmean(Ximg))
    Ximg_std = StandardScaler().fit_transform(Ximg)
    Dp = squareform(pdist(Ximg_std, metric="euclidean"))

    # ---- trait distance (standardised NUE traits, Control scenario) ----
    trait_cols = [f"{t}_Control" for t in NUE_TRAITS if f"{t}_Control" in coh.columns]
    Xtr = coh[trait_cols].to_numpy(dtype=float)
    col_mean = np.nanmean(Xtr, axis=0)
    inds = np.where(np.isnan(Xtr))
    Xtr[inds] = np.take(col_mean, inds[1])
    Xtr_std = StandardScaler().fit_transform(Xtr)
    Dt = squareform(pdist(Xtr_std, metric="euclidean"))

    # ==== 1-2. Mantel + partial Mantel ====
    r_gp, p_gp = mantel(Dg, Dp)
    r_gt, p_gt = mantel(Dg, Dt)
    r_pt, p_pt = mantel(Dp, Dt)
    r_gp_t, p_gp_t = partial_mantel(Dg, Dp, Dt)
    r_gt_p, p_gt_p = partial_mantel(Dg, Dt, Dp)
    print(f"  Mantel genomic<->phenomic: r={r_gp:.4f} p={p_gp:.4f}")
    print(f"  Mantel genomic<->trait:    r={r_gt:.4f} p={p_gt:.4f}")
    print(f"  Mantel phenomic<->trait:   r={r_pt:.4f} p={p_pt:.4f}")
    print(f"  partial Mantel genomic<->phenomic | trait:  r={r_gp_t:.4f} p={p_gp_t:.4f}")
    print(f"  partial Mantel genomic<->trait | phenomic:  r={r_gt_p:.4f} p={p_gt_p:.4f}")

    pd.DataFrame([
        dict(panel=panel, n=n, comparison="genomic~phenomic", r=round(r_gp, 4), p=round(p_gp, 4)),
        dict(panel=panel, n=n, comparison="genomic~trait", r=round(r_gt, 4), p=round(p_gt, 4)),
        dict(panel=panel, n=n, comparison="phenomic~trait", r=round(r_pt, 4), p=round(p_pt, 4)),
        dict(panel=panel, n=n, comparison="genomic~phenomic|trait(partial)", r=round(r_gp_t, 4), p=round(p_gp_t, 4)),
        dict(panel=panel, n=n, comparison="genomic~trait|phenomic(partial)", r=round(r_gt_p, 4), p=round(p_gt_p, 4)),
    ]).to_csv(TAB / f"concordance_mantel_{panel.lower()}.csv", index=False)

    # ==== 3. Procrustes / PROTEST (genomic PCA vs phenomic PCA) ====
    pca_g = pd.read_csv(TAB / f"pca_{panel.lower()}.csv").rename(columns={"sample": "Taxa"})
    pca_g["Taxa"] = pca_g["Taxa"].astype(str).str.strip()
    pca_g = pca_g.set_index("Taxa").reindex(common)
    Xg_pc = pca_g[["PC1", "PC2", "PC3", "PC4"]].to_numpy()
    Xp_pc = PCA(n_components=4, random_state=SEED).fit_transform(Ximg_std)
    disparity, p_prot = protest(Xg_pc, Xp_pc)
    print(f"  Procrustes/PROTEST genomic-PCA vs phenomic-PCA: disparity(M2)={disparity:.4f} p={p_prot:.4f}")
    pd.DataFrame([dict(panel=panel, n=n, disparity_M2=round(disparity, 4), perm_p=round(p_prot, 4),
                        n_perm=N_PERM)]).to_csv(TAB / f"concordance_procrustes_{panel.lower()}.csv", index=False)

    # ==== 4. Supervised classification: phenomic -> admixture cluster ====
    q = pd.read_csv(TAB / f"admixture_{panel.lower()}_Q.csv")
    q["sample"] = q["sample"].astype(str).str.strip()
    q = q.set_index("sample").reindex(common)
    qcols = [c for c in q.columns if c.startswith("Q")]
    y = q[qcols].to_numpy().argmax(axis=1)
    valid = ~np.isnan(q[qcols].to_numpy()).any(axis=1)
    # Unscaled here deliberately: Ximg_std was standardised on the FULL cohort (fine for the
    # Mantel/PCA uses above, which aren't fold-based), but reusing it for cross-validated
    # classification would leak each held-out fold's rows into the scaler fit. The pipeline
    # below refits StandardScaler on the training fold only, inside cross_val_predict.
    Xc, yc = Ximg[valid], y[valid]

    class_counts = pd.Series(yc).value_counts()
    n_splits = min(5, class_counts.min()) if class_counts.min() >= 2 else 2
    if n_splits >= 2 and len(np.unique(yc)) >= 2:
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=SEED)
        clf = make_pipeline(StandardScaler(),
                             RandomForestClassifier(n_estimators=300, random_state=SEED, n_jobs=-1))
        y_pred = cross_val_predict(clf, Xc, yc, cv=skf)
        acc = accuracy_score(yc, y_pred)
        majority = class_counts.max() / class_counts.sum()
        print(f"  RF phenomic->admixture-cluster: CV accuracy={acc:.4f} (majority-class baseline={majority:.4f}, {n_splits}-fold)")
        pd.DataFrame([dict(panel=panel, n=len(yc), n_clusters=len(np.unique(yc)), n_splits=n_splits,
                            cv_accuracy=round(acc, 4), majority_baseline=round(majority, 4))]).to_csv(
            TAB / f"concordance_classification_{panel.lower()}.csv", index=False)
    else:
        print(f"  SKIP classification for {panel}: insufficient class sizes for CV")

    # ==== 5. Feature-type attribution ====
    groups = {}
    for c in img_cols:
        groups.setdefault(feature_group_of(c), []).append(c)
    attr_rows = []
    for gname, cols in groups.items():
        Xg_feat = StandardScaler().fit_transform(np.nan_to_num(coh[cols].to_numpy(dtype=float),
                                                                  nan=np.nanmean(coh[cols].to_numpy(dtype=float))))
        Dp_g = squareform(pdist(Xg_feat, metric="euclidean"))
        r_g, p_g = mantel(Dg, Dp_g, n_perm=N_PERM)
        attr_rows.append(dict(panel=panel, feature_group=gname, n_features=len(cols),
                               r_genomic=round(r_g, 4), p=round(p_g, 4)))
        print(f"  feature group {gname:16s} ({len(cols):3d} features): genomic Mantel r={r_g:.4f} p={p_g:.4f}")
    pd.DataFrame(attr_rows).to_csv(TAB / f"concordance_feature_attribution_{panel.lower()}.csv", index=False)

    # ==== 6. Structure-as-confounder: phenomic~trait vs phenomic~trait|genomic ====
    r_pt_g, p_pt_g = partial_mantel(Dp, Dt, Dg)
    drop_pct = round(100 * (1 - abs(r_pt_g) / abs(r_pt)), 1) if r_pt != 0 else float("nan")
    print(f"  phenomic~trait (raw) r={r_pt:.4f}  |  phenomic~trait|genomic (partial) r={r_pt_g:.4f}  "
          f"-> {drop_pct:+.1f}% change once genomic structure is controlled for")
    pd.DataFrame([dict(panel=panel, r_raw=round(r_pt, 4), p_raw=round(p_pt, 4),
                        r_partial_given_genomic=round(r_pt_g, 4), p_partial=round(p_pt_g, 4),
                        pct_change=drop_pct)]).to_csv(TAB / f"concordance_confounder_{panel.lower()}.csv", index=False)

print("\n-> tables/concordance_{mantel,procrustes,classification,feature_attribution,confounder}_{set}.csv")
