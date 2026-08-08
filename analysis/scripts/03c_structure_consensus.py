#!/usr/bin/env python
"""03c_structure_consensus.py — Stage 2 wrap-up: do PCA, sNMF admixture, and UMAP agree on
structure? Quantifies cross-method concordance (Adjusted Rand Index, Normalised Mutual
Information) between:
  (a) PCA k-means (k=3, from 02_pca_structure.R's first-pass dominant-cluster read)
  (b) sNMF admixture dominant-component assignment (argmax Q, at the elbow K from 03_admixture.R)
  (c) k-means on the 2D UMAP embedding (k = the same elbow K, for a fair comparison with (b))

This is deliberately the SAME concordance toolkit (ARI/NMI) that Pillar B (09_phenomic_concordance)
will later apply between genomic and phenomic cluster labels -- Stage 2 establishes the baseline
question "do independent genomic structure-inference methods even agree with each other" before
Pillar B asks whether phenomic data agrees with genomic structure.

Inputs: tables/pca_{set}.csv (has PC1..PC6 per sample, no cluster label -- recomputed here),
        tables/admixture_{set}_Q.csv, tables/admixture_bestK.csv, tables/umap_{set}.csv
Output: tables/structure_consensus_{set}.csv, tables/structure_consensus_summary.csv
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import TAB

SEED = 42
PANELS = ["Set1", "Set2"]
bestK = pd.read_csv(TAB / "admixture_bestK.csv").set_index("panel")["bestK_elbow"].to_dict()

all_rows = []
for panel in PANELS:
    print(f"\n=== [{panel}] structure consensus ===")
    k = int(bestK[panel])

    pca = pd.read_csv(TAB / f"pca_{panel.lower()}.csv").rename(columns={"sample": "sample"})
    q = pd.read_csv(TAB / f"admixture_{panel.lower()}_Q.csv")
    um = pd.read_csv(TAB / f"umap_{panel.lower()}.csv")

    df = pca.merge(q, on="sample", how="inner").merge(um, on="sample", how="inner")
    n_common = len(df)
    print(f"  {n_common} accessions with PCA + admixture + UMAP all present "
          f"(pca={len(pca)}, admixture={len(q)}, umap={len(um)})")

    # (a) PCA k-means, k=3 (matches 02_pca_structure.R's first-pass dominant-cluster read)
    pca_km = KMeans(n_clusters=3, n_init=25, random_state=SEED).fit(df[["PC1", "PC2", "PC3", "PC4"]])
    lab_pca3 = pca_km.labels_

    # (a2) PCA k-means at the SAME k as admixture, for a fair like-for-like comparison too
    pca_kmK = KMeans(n_clusters=k, n_init=25, random_state=SEED).fit(df[["PC1", "PC2", "PC3", "PC4"]])
    lab_pcaK = pca_kmK.labels_

    # (b) admixture dominant-component assignment (argmax Q) at the elbow K
    qcols = [c for c in df.columns if c.startswith("Q")]
    lab_admix = df[qcols].to_numpy().argmax(axis=1)

    # (c) k-means on the UMAP embedding, same k as admixture
    umap_km = KMeans(n_clusters=k, n_init=25, random_state=SEED).fit(df[["UMAP1", "UMAP2"]])
    lab_umap = umap_km.labels_

    pairs = [
        ("PCA_kmeans_k3", "admixture_argmaxQ", lab_pca3, lab_admix),
        (f"PCA_kmeans_k{k}", "admixture_argmaxQ", lab_pcaK, lab_admix),
        (f"UMAP_kmeans_k{k}", "admixture_argmaxQ", lab_umap, lab_admix),
        (f"PCA_kmeans_k{k}", f"UMAP_kmeans_k{k}", lab_pcaK, lab_umap),
    ]
    for name_a, name_b, la, lb in pairs:
        ari = adjusted_rand_score(la, lb)
        nmi = normalized_mutual_info_score(la, lb)
        print(f"  {name_a:20s} vs {name_b:20s}  ARI={ari:6.3f}  NMI={nmi:6.3f}")
        all_rows.append(dict(panel=panel, n=n_common, k_used=k,
                              method_a=name_a, method_b=name_b, ARI=round(ari, 4), NMI=round(nmi, 4)))

    pd.DataFrame([r for r in all_rows if r["panel"] == panel]).to_csv(
        TAB / f"structure_consensus_{panel.lower()}.csv", index=False)

out = pd.DataFrame(all_rows)
out.to_csv(TAB / "structure_consensus_summary.csv", index=False)
print(f"\n-> tables/structure_consensus_{{set1,set2}}.csv + structure_consensus_summary.csv")
