#!/usr/bin/env python
"""32_overlap_null.py — [stage 32] a chance expectation for the
"66% of haplotype-scan regions overlap pcadapt outliers" convergence claim, plus the pcadapt
saturation and array-coverage bookkeeping.

Null model: per replicate, the haplotype-scan outlier SNP positions are CIRCULARLY SHIFTED
within each chromosome by an independent uniform offset (preserving each chromosome's outlier
count and clustering structure), regions re-merged with the same 200-kb rule, and the pcadapt-
overlap fraction (+/-100 kb, as in stage 17) recomputed. 999 rotations -> empirical P and
enrichment ratio (observed / null mean). An exact-SNP overlap statistic is reported alongside:
the fraction of iHS/XP-EHH outlier SNPs that are themselves pcadapt outliers vs the pcadapt
outlier density (binomial expectation).

Bookkeeping: count of Set 2 pcadapt P-values saturated at the numerical floor, and per-
chromosome QC SNP counts for both panels (the chr4/chr7 'deserts' in Fig 9b are array-design
coverage, shown directly).

Outputs: overlap_null_summary.csv, pcadapt_saturation.csv, snp_density_by_chrom.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import TAB, INTERIM

RNG = np.random.default_rng(42)
N_ROT = 999
MERGE_GAP = 200_000
FLANK = 100_000


def merge_regions(chroms, poss):
    df = pd.DataFrame({"chrom": chroms, "pos": poss}).sort_values(["chrom", "pos"])
    regions = []
    for ch, sub in df.groupby("chrom"):
        p = sub["pos"].to_numpy()
        start = p[0]; prev = p[0]; n = 1
        for x in p[1:]:
            if x - prev > MERGE_GAP:
                regions.append((ch, start, prev, n))
                start, n = x, 0
            prev = x; n += 1
        regions.append((ch, start, prev, n))
    return pd.DataFrame(regions, columns=["chrom", "start", "end", "n_snps"])


def overlap_fraction(regions, pc_by_chrom):
    hits = 0
    for _, r in regions.iterrows():
        pc = pc_by_chrom.get(r["chrom"])
        if pc is not None and np.any((pc >= r["start"] - FLANK) & (pc <= r["end"] + FLANK)):
            hits += 1
    return hits / len(regions), hits


def main():
    ihs = pd.read_csv(TAB / "ihs_set1.csv")
    xp = pd.read_csv(TAB / "xpehh_set1.csv")
    hap = pd.concat([
        ihs.loc[ihs["outlier"], ["chrom", "pos"]],
        xp.loc[xp["outlier"], ["chrom", "pos"]]]).astype({"chrom": str})
    pc = pd.read_csv(TAB / "pcadapt_outliers_set1.csv")
    pc["chrom"] = pc["chrom"].astype(str)
    pc_out = pc[pc["outlier"]]
    pc_by_chrom = {ch: sub["pos"].to_numpy() for ch, sub in pc_out.groupby("chrom")}
    chrom_len = pc.groupby("chrom")["pos"].max().to_dict()

    obs_regions = merge_regions(hap["chrom"], hap["pos"])
    obs_frac, obs_hits = overlap_fraction(obs_regions, pc_by_chrom)
    print(f"[32] observed: {len(obs_regions)} regions, {obs_hits} overlap pcadapt "
          f"({100*obs_frac:.1f}%)")

    null_fracs = []
    for _ in range(N_ROT):
        rot_ch, rot_pos = [], []
        for ch, sub in hap.groupby("chrom"):
            L = chrom_len.get(ch)
            if L is None:
                continue
            off = RNG.integers(1, L)
            rot_ch.extend([ch] * len(sub))
            rot_pos.extend(((sub["pos"].to_numpy() + off) % L).tolist())
        regs = merge_regions(np.array(rot_ch), np.array(rot_pos))
        f, _ = overlap_fraction(regs, pc_by_chrom)
        null_fracs.append(f)
    null_fracs = np.array(null_fracs)
    p_emp = float((np.sum(null_fracs >= obs_frac) + 1) / (N_ROT + 1))
    enrich = float(obs_frac / null_fracs.mean())

    # exact-SNP overlap: haplotype outlier SNPs that are ALSO pcadapt outliers
    pc_set = set(zip(pc_out["chrom"], pc_out["pos"]))
    hap_pairs = list(zip(hap["chrom"], hap["pos"]))
    exact = sum(1 for t in hap_pairs if t in pc_set)
    dens = len(pc_out) / len(pc)
    exact_enrich = (exact / len(hap_pairs)) / dens

    pd.DataFrame([{
        "n_regions": len(obs_regions), "n_overlap": obs_hits,
        "observed_fraction": obs_frac, "null_mean_fraction": float(null_fracs.mean()),
        "null_p97_5": float(np.percentile(null_fracs, 97.5)),
        "empirical_p": p_emp, "enrichment": enrich,
        "n_singleton_regions": int((obs_regions["n_snps"] == 1).sum()),
        "exact_snp_overlap": exact, "n_hap_outlier_snps": len(hap_pairs),
        "exact_overlap_fraction": exact / len(hap_pairs),
        "pcadapt_outlier_density": dens, "exact_snp_enrichment": exact_enrich,
        "n_rotations": N_ROT}]).to_csv(TAB / "overlap_null_summary.csv", index=False)
    print(f"[32] rotation null: mean {100*null_fracs.mean():.1f}%, "
          f"P = {p_emp:.4f}, enrichment = {enrich:.2f}x")
    print(f"[32] exact-SNP: {exact}/{len(hap_pairs)} = "
          f"{100*exact/len(hap_pairs):.1f}% vs density {100*dens:.2f}% "
          f"-> {exact_enrich:.1f}x")

    # pcadapt saturation + per-chromosome SNP density
    sat_rows = []
    for panel in ["set1", "set2"]:
        d = pd.read_csv(TAB / f"pcadapt_outliers_{panel}.csv")
        pmin = d["pvalue"][d["pvalue"] > 0].min()
        n_floor = int((d["pvalue"] <= pmin).sum())
        n_zero = int((d["pvalue"] == 0).sum())
        sat_rows.append({"panel": panel.title(), "min_positive_p": pmin,
                         "n_at_floor": n_floor, "n_exact_zero": n_zero,
                         "n_snps": len(d)})
    pd.DataFrame(sat_rows).to_csv(TAB / "pcadapt_saturation.csv", index=False)

    dens_rows = []
    for panel in ["set1", "set2"]:
        bim = pd.read_csv(INTERIM / f"{panel.title()}_qc.bim", sep=r"\s+", header=None,
                          names=["chrom", "snp", "cm", "pos", "a1", "a2"])
        g = bim.groupby("chrom").agg(n_snps=("snp", "size"), span_mb=("pos", "max"))
        g["snps_per_mb"] = g["n_snps"] / (g["span_mb"] / 1e6)
        g = g.reset_index(); g.insert(0, "panel", panel.title())
        dens_rows.append(g)
    pd.concat(dens_rows).to_csv(TAB / "snp_density_by_chrom.csv", index=False)
    print("[32] done")


if __name__ == "__main__":
    main()
