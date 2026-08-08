#!/usr/bin/env python
"""06_diversity.py — classic diversity stats (MAF spectrum, PIC, He/Ho/inbreeding F via PLINK)
+ genome-wide sliding-window nucleotide diversity (pi, Watterson's theta, Tajima's D) via
scikit-allel, computed on the QC (full, NOT LD-pruned) SNP set -- pruning removes physically
clustered SNPs and would bias window-based density estimates.

Ne (effective population size, LD-based) is explicitly NOT computed here: a correct LD-based Ne
estimator (e.g. SNeP-style) is non-trivial to reimplement correctly from scratch, and an
incorrectly implemented number would be worse than an honest gap. Deferred; note this in the
manuscript's diversity section rather than reporting an unverified figure.

Outputs (analysis/results/tables/):
  maf_spectrum_{set}.csv       — MAF bin, n_snps
  diversity_summary_{set}.csv  — one row: n, n_snps, mean_PIC, mean_He, mean_Ho, mean_F,
                                  genome_pi, genome_theta_w, genome_tajima_d
  diversity_windows_{set}.csv  — per-window (chrom, window_start, window_end, pi, theta_w, tajima_d)
"""
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import allel

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PLINK, INTERIM, TAB

WINDOW_BP = 1_000_000
PANELS = ["Set1", "Set2"]


def run_plink(*args):
    cmd = [str(PLINK), *[str(a) for a in args], "--chr-set", "12", "no-xy",
           "--allow-extra-chr", "--silent"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("PLINK ERROR:\n", r.stderr[-1500:])
        raise SystemExit(1)
    return r


def maf_pic_stats(frq_path):
    frq = pd.read_csv(frq_path, sep=r"\s+")
    p = frq["MAF"].to_numpy(dtype=float)
    q = 1 - p
    pic = 1 - p**2 - q**2 - 2 * p**2 * q**2
    bins = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
    labels = ["0-0.05", "0.05-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5"]
    cats = pd.cut(p, bins=bins, labels=labels, include_lowest=True)
    spectrum = cats.value_counts().reindex(labels).reset_index()
    spectrum.columns = ["maf_bin", "n_snps"]
    return spectrum, float(np.nanmean(pic)), len(frq)


def het_stats(het_path):
    het = pd.read_csv(het_path, sep=r"\s+")
    n = het["N(NM)"].to_numpy(dtype=float)
    ohom = het["O(HOM)"].to_numpy(dtype=float)
    ehom = het["E(HOM)"].to_numpy(dtype=float)
    ho = 1 - ohom / n
    he = 1 - ehom / n
    f = het["F"].to_numpy(dtype=float)
    return float(np.nanmean(ho)), float(np.nanmean(he)), float(np.nanmean(f))


def windowed_diversity(vcf_path):
    callset = allel.read_vcf(str(vcf_path), fields=["variants/CHROM", "variants/POS", "calldata/GT"])
    chrom = callset["variants/CHROM"]
    pos = callset["variants/POS"]
    gt = allel.GenotypeArray(callset["calldata/GT"])
    ac = gt.count_alleles()

    rows = []
    chrom_pi_rows = []   # per-chromosome genome-wide pi/theta, to combine into a genome-wide figure
                          # WITHOUT concatenating positions across chromosomes (they are not
                          # monotonically increasing across chromosome boundaries -- each
                          # chromosome restarts near position 1, which crashed the naive
                          # genome-concatenated allel.sequence_diversity() call; fixed by computing
                          # per-chromosome then combining with a site-count-weighted mean, the
                          # standard way to report a genome-wide average of a per-site statistic)
    for c in sorted(set(chrom), key=lambda x: (len(x), x)):
        mask = chrom == c
        p = pos[mask]
        a = ac[mask]
        if len(p) < 10:
            continue
        order = np.argsort(p)
        p, a = p[order], a[order]
        pi, windows, n_bases, counts = allel.windowed_diversity(p, a, size=WINDOW_BP, start=p.min(), stop=p.max())
        theta, _, _, _ = allel.windowed_watterson_theta(p, a, size=WINDOW_BP, start=p.min(), stop=p.max())
        tajd, _, _ = allel.windowed_tajima_d(p, a, size=WINDOW_BP, start=p.min(), stop=p.max(), min_sites=3)
        for i in range(len(windows)):
            rows.append(dict(chrom=c, window_start=int(windows[i, 0]), window_end=int(windows[i, 1]),
                              n_snps=int(counts[i]), pi=pi[i], theta_w=theta[i], tajima_d=tajd[i]))

        chrom_pi = float(allel.sequence_diversity(p, a))
        chrom_theta = float(allel.watterson_theta(p, a))
        chrom_pi_rows.append(dict(chrom=c, n_snps=len(p), pi=chrom_pi, theta_w=chrom_theta))

    win_df = pd.DataFrame(rows)
    cpr = pd.DataFrame(chrom_pi_rows)
    # site-count-weighted genome-wide average across chromosomes (not a naive concatenation)
    genome_pi = float(np.average(cpr["pi"], weights=cpr["n_snps"])) if len(cpr) else float("nan")
    genome_theta = float(np.average(cpr["theta_w"], weights=cpr["n_snps"])) if len(cpr) else float("nan")
    genome_tajd = float(win_df["tajima_d"].mean()) if len(win_df) else float("nan")
    return win_df, genome_pi, genome_theta, genome_tajd


rows_summary = []
for panel in PANELS:
    print(f"\n=== [{panel}] diversity ===")
    stem = INTERIM / f"{panel}_qc"
    out_stub = INTERIM / "lea" / f"{panel}_diversity"
    out_stub.parent.mkdir(parents=True, exist_ok=True)

    # ---- PLINK: allele frequencies (MAF spectrum + PIC) ----
    run_plink("--bfile", stem, "--freq", "--out", out_stub)
    spectrum, mean_pic, n_snps = maf_pic_stats(out_stub.with_suffix(".frq"))
    spectrum.insert(0, "panel", panel)
    spectrum.to_csv(TAB / f"maf_spectrum_{panel.lower()}.csv", index=False)
    print(f"  MAF spectrum + PIC: {n_snps} SNPs, mean PIC = {mean_pic:.4f}")

    # ---- PLINK: heterozygosity / inbreeding F ----
    run_plink("--bfile", stem, "--het", "--out", out_stub)
    het_df = pd.read_csv(out_stub.with_suffix(".het"), sep=r"\s+")
    n_samples = len(het_df)
    mean_ho, mean_he, mean_f = het_stats(out_stub.with_suffix(".het"))
    print(f"  n={n_samples}  Ho={mean_ho:.4f}  He={mean_he:.4f}  F={mean_f:.4f}")

    # ---- scikit-allel: windowed pi / theta_w / Tajima's D ----
    vcf_out = out_stub.with_name(out_stub.name + "_vcf")
    run_plink("--bfile", stem, "--recode", "vcf", "--out", vcf_out)
    win_df, genome_pi, genome_theta, genome_tajd = windowed_diversity(vcf_out.with_suffix(".vcf"))
    win_df.insert(0, "panel", panel)
    win_df.to_csv(TAB / f"diversity_windows_{panel.lower()}.csv", index=False)
    print(f"  genome-wide pi={genome_pi:.5f}  theta_w={genome_theta:.5f}  mean windowed TajD={genome_tajd:.4f}")
    print(f"  {len(win_df)} windows ({WINDOW_BP/1e6:.0f} Mb) across {win_df['chrom'].nunique()} chromosomes")

    rows_summary.append(dict(panel=panel, n=n_samples, n_snps=n_snps, mean_PIC=round(mean_pic, 4),
                              mean_Ho=round(mean_ho, 4), mean_He=round(mean_he, 4), mean_F=round(mean_f, 4),
                              genome_pi=round(genome_pi, 6), genome_theta_w=round(genome_theta, 6),
                              genome_tajima_d=round(genome_tajd, 4),
                              Ne_LD_based="not computed (deferred, see script header)"))

pd.DataFrame(rows_summary).to_csv(TAB / "diversity_summary.csv", index=False)
print("\n-> tables/maf_spectrum_{set}.csv + diversity_windows_{set}.csv + diversity_summary.csv")
