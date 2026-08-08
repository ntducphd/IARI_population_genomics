#!/usr/bin/env python
"""16_hwe_selfing.py — [stage 16] Hardy-Weinberg departures, Fis distribution, and a
mating-system (selfing-rate) estimate per panel.

Motivation: rice is predominantly selfing, so every
HWE-derived quantity in the paper (He, F, Gst) operates far from panmixia. Rather than leaving F
as an unused table entry, this stage makes the mating system explicit: the equilibrium
inbreeding-selfing relation F = s/(2-s) inverts to s = 2F/(1+F) (Wright's equilibrium under
partial self-fertilisation), giving an interpretable selfing-rate estimate, and the HWE exact-test
spectrum documents how completely heterozygote deficit dominates the panels.

Caveats stated in the manuscript, not hidden: s-hat assumes inbreeding equilibrium and no
selection/assortative structure, and F here averages over accessions maintained as inbred lines —
it estimates the realised, breeding-programme mating system, not the outcrossing biology of wild
rice.

Outputs (analysis/results/tables/):
  hwe_summary.csv         — per panel: n_snps_tested, share P<1e-6 (het deficit), share het excess,
                             mean_Fis, median_Fis, selfing_rate_hat (from mean and median F)
  fis_per_accession_{set}.csv — sample, cluster, F (PLINK --het), s_hat
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PLINK, INTERIM, TAB

PANELS = ["Set1", "Set2"]


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


def main():
    rows = []
    for panel in PANELS:
        bfile = INTERIM / f"{panel}_qc"
        out = INTERIM / f"{panel}_hwe"
        run_plink("--bfile", bfile, "--hardy", "--het", "--out", out)

        hwe = pd.read_csv(f"{out}.hwe", sep=r"\s+")
        # PLINK 1.9 labels the all-samples test "ALL" with a phenotype and "ALL(NP)" without one
        hwe = hwe[hwe["TEST"].str.startswith("ALL")] if "TEST" in hwe.columns else hwe
        # direction of departure: observed vs expected heterozygosity per SNP
        het_deficit = hwe["O(HET)"] < hwe["E(HET)"]
        sig = hwe["P"] < 1e-6

        het = pd.read_csv(f"{out}.het", sep=r"\s+")
        f = het["F"].to_numpy(dtype=float)
        s_mean = 2 * np.nanmean(f) / (1 + np.nanmean(f))
        s_median = 2 * np.nanmedian(f) / (1 + np.nanmedian(f))

        clu = cluster_map(panel)
        acc = pd.DataFrame({"sample": het["IID"], "F": f})
        acc["cluster"] = acc["sample"].map(clu)
        acc["s_hat"] = 2 * acc["F"] / (1 + acc["F"])
        acc.to_csv(TAB / f"fis_per_accession_{panel.lower()}.csv", index=False)

        rows.append({
            "panel": panel, "n_snps_tested": len(hwe),
            "share_sig_het_deficit": float((sig & het_deficit).mean()),
            "share_sig_het_excess": float((sig & ~het_deficit).mean()),
            "mean_F": float(np.nanmean(f)), "median_F": float(np.nanmedian(f)),
            "selfing_rate_from_mean_F": s_mean, "selfing_rate_from_median_F": s_median})
        print(f"[16] {panel}: mean F = {np.nanmean(f):.3f} -> s_hat = {s_mean:.3f}; "
              f"{(sig & het_deficit).mean():.1%} of SNPs in significant het deficit")

    pd.DataFrame(rows).to_csv(TAB / "hwe_summary.csv", index=False)
    print("[16] done -> hwe_summary.csv")


if __name__ == "__main__":
    main()
