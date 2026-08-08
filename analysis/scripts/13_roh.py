#!/usr/bin/env python
"""13_roh.py — [stage 13] runs of homozygosity (ROH) and genomic inbreeding F_ROH per panel.

Motivation: ROH is the standard genomic-inbreeding
measure in human and livestock genetics; its length spectrum separates recent from ancient
inbreeding. In a predominantly selfing crop the expectation itself is the teaching point: repeated
selfing drives most of the genome into (long) ROH, so F_ROH should approach 1 and the informative
signal is the length-class spectrum and the residual heterozygosity islands — the opposite regime
from the F_ROH ~ 0.01-0.3 typical of outbreeding humans/livestock.

Method: PLINK 1.9 --homozyg on the QC (full, NOT LD-pruned) SNP set — pruning would delete the
physically clustered SNPs ROH detection depends on. Parameters are set per panel because marker
density differs ~16x (Set1 WGS-derived ~1.35 SNP/kb; Set2 50K array ~1 SNP/12 kb); Set2 uses
relaxed, livestock-50K-style settings. Genome length for F_ROH = sum over chromosomes of the max
SNP position in the panel's own .bim (a slight underestimate of true chromosome ends, applied
identically to both panels).

Outputs (analysis/results/tables/):
  roh_indiv_{set}.csv          — per accession: cluster, n_roh, total_kb, mean_kb, f_roh
  roh_length_classes_{set}.csv — per length class (<0.5, 0.5-1, 1-2, 2-4, 4-8, >8 Mb):
                                  n_segments, total_mb, share_of_roh
  roh_cluster_summary_{set}.csv— per admixture cluster: n, mean_f_roh, sd_f_roh, mean_n_roh
  roh_summary.csv              — one row per panel: n, mean/median F_ROH, mean total ROH Mb,
                                  genome_mb used, parameter string (provenance)
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PLINK, INTERIM, TAB

PANELS = {
    # window/segment parameters chosen per marker density; recorded verbatim in roh_summary.csv
    "Set1": ["--homozyg-window-snp", "50", "--homozyg-window-het", "1",
             "--homozyg-window-missing", "5", "--homozyg-snp", "100",
             "--homozyg-kb", "300", "--homozyg-density", "50", "--homozyg-gap", "1000"],
    "Set2": ["--homozyg-window-snp", "20", "--homozyg-window-het", "1",
             "--homozyg-window-missing", "3", "--homozyg-snp", "20",
             "--homozyg-kb", "1000", "--homozyg-density", "200", "--homozyg-gap", "1000"],
}
LEN_BINS = [0, 500, 1000, 2000, 4000, 8000, np.inf]              # kb
LEN_LABELS = ["<0.5Mb", "0.5-1Mb", "1-2Mb", "2-4Mb", "4-8Mb", ">8Mb"]


def run_plink(*args):
    cmd = [str(PLINK), *[str(a) for a in args], "--chr-set", "12", "no-xy",
           "--allow-extra-chr", "--silent"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("PLINK ERROR:\n", r.stderr[-1500:])
        raise SystemExit(1)


def cluster_map(panel):
    q = pd.read_csv(TAB / f"admixture_{panel.lower()}_Q.csv")
    qcols = [c for c in q.columns if c.startswith("Q")]
    return dict(zip(q["sample"], "C" + q[qcols].to_numpy().argmax(axis=1).astype(str)))


def genome_kb(bim_path):
    bim = pd.read_csv(bim_path, sep=r"\s+", header=None,
                      names=["chrom", "snp", "cm", "pos", "a1", "a2"])
    return float(bim.groupby("chrom")["pos"].max().sum()) / 1000.0


def main():
    summary_rows = []
    for panel, params in PANELS.items():
        bfile = INTERIM / f"{panel}_qc"
        out = INTERIM / f"{panel}_roh"
        run_plink("--bfile", bfile, "--homozyg", *params, "--out", out)

        hom = pd.read_csv(f"{out}.hom", sep=r"\s+")
        indiv = pd.read_csv(f"{out}.hom.indiv", sep=r"\s+")
        gkb = genome_kb(f"{bfile}.bim")
        clu = cluster_map(panel)

        indiv = indiv.rename(columns={"IID": "sample", "NSEG": "n_roh", "KB": "total_kb"})
        indiv["mean_kb"] = indiv["total_kb"] / indiv["n_roh"].replace(0, np.nan)
        indiv["f_roh"] = indiv["total_kb"] / gkb
        indiv["cluster"] = indiv["sample"].map(clu)
        indiv[["sample", "cluster", "n_roh", "total_kb", "mean_kb", "f_roh"]].to_csv(
            TAB / f"roh_indiv_{panel.lower()}.csv", index=False)

        cats = pd.cut(hom["KB"], bins=LEN_BINS, labels=LEN_LABELS, include_lowest=True)
        lc = hom.groupby(cats, observed=False)["KB"].agg(n_segments="size", total_kb="sum")
        lc["total_mb"] = lc["total_kb"] / 1000.0
        lc["share_of_roh"] = lc["total_kb"] / lc["total_kb"].sum()
        lc.drop(columns="total_kb").reset_index(names="length_class").to_csv(
            TAB / f"roh_length_classes_{panel.lower()}.csv", index=False)

        cs = indiv.groupby("cluster", dropna=False).agg(
            n=("sample", "size"), mean_f_roh=("f_roh", "mean"),
            sd_f_roh=("f_roh", "std"), mean_n_roh=("n_roh", "mean")).reset_index()
        cs.to_csv(TAB / f"roh_cluster_summary_{panel.lower()}.csv", index=False)

        summary_rows.append({
            "panel": panel, "n": len(indiv), "n_segments": len(hom),
            "mean_f_roh": indiv["f_roh"].mean(), "median_f_roh": indiv["f_roh"].median(),
            "mean_total_roh_mb": indiv["total_kb"].mean() / 1000.0,
            "genome_mb": gkb / 1000.0, "plink_params": " ".join(params)})
        print(f"[13] {panel}: {len(hom)} segments, "
              f"mean F_ROH = {indiv['f_roh'].mean():.3f}, "
              f"median = {indiv['f_roh'].median():.3f}")

    pd.DataFrame(summary_rows).to_csv(TAB / "roh_summary.csv", index=False)
    print("[13] done -> roh_summary.csv")


if __name__ == "__main__":
    main()
