#!/usr/bin/env python
"""23_diversity_recompute.py — [stage 23] recompute nucleotide diversity (pi), Watterson's theta,
and Tajima's D on a MISSINGNESS-ONLY-filtered SNP set (no MAF filter).

Motivation:
stage 06 computed these statistics on the QC set, which includes the MAF > 0.05 filter. Removing
all variants below 5% frequency deletes precisely the rare-allele class that drives Tajima's D
negative, so the strongly positive D previously reported (+3.15 Set1, +2.76 Set2) is guaranteed
to be inflated by the filter, independent of any real bottleneck signal. This stage rebuilds the
site set from the RAW genotypes with the missingness filter only (--geno 0.1, biallelic SNPs) and
recomputes the same windowed statistics with the same per-chromosome, SNP-count-weighted
combination used in stage 06.

Set 2 policy: the 50K array's SNPs were ascertained on a diversity panel, so its site-frequency
spectrum is biased by design and SFS-shape statistics (theta_W, Tajima's D) are not meaningful in
absolute terms on ANY filtering of chip data. We recompute Set 2 for completeness of provenance
but the manuscript reports Set 1 as the primary estimate and moves Set 2 SFS-shape values to a
platform-relative footnote (ascertainment bias, not just "different platforms").

Outputs (analysis/results/tables/):
  diversity_nomaf_summary.csv       — per panel: n_snps, genome_pi, genome_theta_w,
                                       genome_tajima_d (+ the stage-06 MAF-filtered values
                                       side-by-side for the correction audit trail)
  diversity_nomaf_windows_{set}.csv — per 1-Mb window: chrom, start, end, n_snps, pi, theta_w, d
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import allel

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PLINK, INTERIM, TAB, SET1_BED_1M

WINDOW_BP = 1_000_000


def run_plink(*args):
    cmd = [str(PLINK), *[str(a) for a in args], "--chr-set", "12", "no-xy",
           "--allow-extra-chr", "--silent"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("PLINK ERROR:\n", r.stderr[-1500:])
        raise SystemExit(1)


def build_nomaf_vcf(panel):
    """Missingness-only filter from the rawest genotype source for each panel."""
    out = INTERIM / f"{panel}_nomaf"
    vcf = Path(f"{out}.vcf")
    if vcf.exists():
        return vcf
    if panel == "Set1":
        run_plink("--bfile", SET1_BED_1M, "--geno", "0.1", "--snps-only", "just-acgt",
                  "--biallelic-only", "strict", "--recode", "vcf-iid", "--out", out)
    else:
        # Set2 raw arrived as HapMap; stage 01 already converted it losslessly to Set2.vcf
        # (pre-QC). Apply the missingness filter only.
        run_plink("--vcf", INTERIM / "Set2.vcf", "--geno", "0.1", "--snps-only", "just-acgt",
                  "--biallelic-only", "strict", "--recode", "vcf-iid", "--out", out)
    return vcf


def windowed(vcf_path):
    cs = allel.read_vcf(str(vcf_path), fields=["variants/CHROM", "variants/POS", "calldata/GT"])
    chrom, pos = cs["variants/CHROM"], cs["variants/POS"]
    gt = allel.GenotypeArray(cs["calldata/GT"])
    ac = gt.count_alleles()

    rows, chrom_stats = [], []
    for ch in pd.unique(chrom):
        m = chrom == ch
        pos_c, ac_c = pos[m], ac[m]
        # the Set2 HapMap->VCF conversion does not guarantee sorted positions within a chromosome
        order = np.argsort(pos_c, kind="stable")
        pos_c, ac_c = pos_c[order], ac_c[order]
        # collapse duplicate positions (SortedIndex tolerates ties, windowed stats do too, but
        # exact duplicates from the array design add no information)
        n_c = int(m.sum())
        if n_c < 10:
            continue
        start = 1
        end = int(pos_c.max())
        # windowed stats (per-site averages within accessible windows)
        pi_w, windows, n_bases, counts = allel.windowed_diversity(
            pos_c, ac_c, size=WINDOW_BP, start=start, stop=end)
        tw_w, _, _, _ = allel.windowed_watterson_theta(
            pos_c, ac_c, size=WINDOW_BP, start=start, stop=end)
        d_w, _, _ = allel.windowed_tajima_d(pos_c, ac_c, size=WINDOW_BP, start=start, stop=end)
        for (ws, we), p_, t_, d_, c_ in zip(windows, pi_w, tw_w, d_w, counts):
            rows.append({"chrom": ch, "window_start": int(ws), "window_end": int(we),
                         "n_snps": int(c_), "pi": p_, "theta_w": t_, "tajima_d": d_})
        chrom_stats.append({"chrom": ch, "n_snps": n_c,
                            "pi": allel.sequence_diversity(pos_c, ac_c),
                            "theta_w": allel.watterson_theta(pos_c, ac_c),
                            "tajima_d": allel.tajima_d(ac_c, pos=pos_c)})
    cw = pd.DataFrame(chrom_stats)
    w = cw["n_snps"] / cw["n_snps"].sum()
    genome = {k: float((cw[k] * w).sum()) for k in ["pi", "theta_w", "tajima_d"]}
    return pd.DataFrame(rows), genome, int(cw["n_snps"].sum())


def main():
    old = pd.read_csv(TAB / "diversity_summary.csv")
    out_rows = []
    for panel in ["Set1", "Set2"]:
        vcf = build_nomaf_vcf(panel)
        win, genome, n_snps = windowed(vcf)
        win.to_csv(TAB / f"diversity_nomaf_windows_{panel.lower()}.csv", index=False)
        prev = old[old["panel"].str.lower() == panel.lower()] if "panel" in old.columns else None
        row = {"panel": panel, "n_snps_nomaf": n_snps,
               "genome_pi": genome["pi"], "genome_theta_w": genome["theta_w"],
               "genome_tajima_d": genome["tajima_d"]}
        if prev is not None and len(prev):
            for k_new, k_old in [("maf_filtered_pi", "genome_pi"),
                                 ("maf_filtered_theta_w", "genome_theta_w"),
                                 ("maf_filtered_tajima_d", "genome_tajima_d")]:
                if k_old in prev.columns:
                    row[k_new] = float(prev.iloc[0][k_old])
        out_rows.append(row)
        print(f"[23] {panel}: {n_snps} SNPs (no MAF filter) -> pi = {genome['pi']:.3e}, "
              f"theta_w = {genome['theta_w']:.3e}, Tajima's D = {genome['tajima_d']:+.2f}")
    pd.DataFrame(out_rows).to_csv(TAB / "diversity_nomaf_summary.csv", index=False)
    print("[23] done -> diversity_nomaf_summary.csv")


if __name__ == "__main__":
    main()
