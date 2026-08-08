#!/usr/bin/env python
"""03b_umap_structure.py — non-linear structure embedding (UMAP), the third independent
triangulation method for Stage 2 (alongside PCA in 02_pca_structure.R and sNMF admixture in
03_admixture.R). DAPC (adegenet) was planned as a fourth method but is DROPPED: adegenet's
read.PLINK crashed R deterministically on this machine when its function signature was even
inspected (environment audit -- likely a broken compiled
dependency given Rtools 4.4 is not installed here). PCA + sNMF + UMAP is still a defensible
triangulation (linear ordination + model-based clustering + non-linear embedding).

Pipeline: Set{1,2}_pruned.{bed,bim,fam} -> PLINK --recode A (additive dosage 0/1/2)
  -> mean-impute missing -> UMAP(n_components=2) -> umap_{set}.csv (sample, UMAP1, UMAP2)
"""
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import umap

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PLINK, INTERIM, TAB

SEED = 42
PANELS = ["Set1", "Set2"]

for panel in PANELS:
    print(f"\n=== [{panel}] UMAP embedding ===")
    stem = INTERIM / f"{panel}_pruned"
    raw_out = INTERIM / "lea" / f"{panel}_umap_raw"
    raw_out.parent.mkdir(parents=True, exist_ok=True)

    cmd = [str(PLINK), "--bfile", str(stem), "--recode", "A",
           "--chr-set", "12", "no-xy", "--allow-extra-chr", "--silent",
           "--out", str(raw_out)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    raw_file = raw_out.with_suffix(".raw")
    if r.returncode != 0 or not raw_file.exists():
        print("PLINK ERROR:\n", r.stderr[-1500:])
        raise SystemExit(1)

    df = pd.read_csv(raw_file, sep=r"\s+")
    meta_cols = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]
    snp_cols = [c for c in df.columns if c not in meta_cols]
    X = df[snp_cols].to_numpy(dtype=float)
    col_mean = np.nanmean(X, axis=0)
    inds = np.where(np.isnan(X))
    X[inds] = np.take(col_mean, inds[1])
    print(f"  {X.shape[0]} accessions x {X.shape[1]} pruned SNPs (missing mean-imputed)")

    n_neighbors = min(15, max(2, X.shape[0] - 1))
    reducer = umap.UMAP(n_components=2, n_neighbors=n_neighbors, min_dist=0.3,
                         random_state=SEED)
    emb = reducer.fit_transform(X)

    out = pd.DataFrame({"sample": df["IID"], "UMAP1": emb[:, 0], "UMAP2": emb[:, 1]})
    out_path = TAB / f"umap_{panel.lower()}.csv"
    out.to_csv(out_path, index=False)
    print(f"  -> {out_path}")

print("\n-> tables/umap_{set1,set2}.csv")
