#!/usr/bin/env python
"""31_demography_repair.py — [stage 31] genome-history layer hardening: make the layer
structure-aware and numerically clean.

(a) Stairway output cleaned: re-parse the final summary, drop numeric-underflow rows
    (year < 0.1) and consecutive duplicates, rewrite tables/stairway_ne_set1.csv (SuppFig S4 is
    refit by 28_figures_flagship.py from the cleaned table).
(b) WITHIN-CLUSTER companion statistics (the structure-controlled numbers the reviews demand):
    for the largest Set 1 admixture cluster — windowed Tajima's D on the missingness-only site
    set restricted to cluster members, and an LD-based Ne from a within-cluster r^2 curve
    (PLINK --keep, same thinning/binning as stage 07, Sved estimator as stage 14 with the
    cluster's own F). If D within the largest cluster is far below the pooled +1.14, the pooled
    value is (as the reviews argue) substantially a structure artefact — reported either way.
(c) Ne INTERVAL reporting: recent Ne recomputed under both inbreeding choices (PLINK F vs
    F_ROH) and the 3-5 cM/Mb map band -> tables/ne_interval_summary.csv, with the effective
    epoch t = 1/(2c*) stated per estimate.

Outputs: stairway_ne_set1.csv (cleaned, overwrites), within_cluster_history_set1.csv,
ne_interval_summary.csv
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import allel

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PLINK, INTERIM, TAB

CM_GRID = (3.0, 4.0, 5.0)
MIN_BIN_BP = 25_000
ALPHA = 1.0


def run_plink(*args):
    cmd = [str(PLINK), *[str(a) for a in args], "--chr-set", "12", "no-xy",
           "--allow-extra-chr", "--silent"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("PLINK ERROR:\n", r.stderr[-1500:])
        raise SystemExit(1)


# ---------- (a) clean the Stairway summary ----------
def clean_stairway():
    final = list((INTERIM / "stairway/set1_out").glob("*final.summary"))
    if not final:
        print("[31a] no Stairway summary found; keeping existing table")
        return
    fs = pd.read_csv(final[0], sep=r"\s+")
    cols = {c.lower(): c for c in fs.columns}
    year_c = cols.get("year")
    # the summary carries BOTH theta_per_site_* and Ne_* columns — select the Ne_ family
    # explicitly (the first-"median"-match heuristic picked theta and yielded Ne = 0)
    ne_cols = [c for c in fs.columns if c.lower().startswith("ne")]
    med_c = [c for c in ne_cols if "median" in c.lower()][0]
    lo_c = [c for c in ne_cols if "2.5" in c][0]
    hi_c = [c for c in ne_cols if "97.5" in c][0]
    out = fs[[year_c, med_c, lo_c, hi_c]].copy()
    out.columns = ["year", "ne_median", "ne_lo95", "ne_hi95"]
    out = out[np.isfinite(out["year"]) & (out["year"] >= 0.1)]
    out = out.loc[(out[["ne_median", "ne_lo95", "ne_hi95"]].diff().abs().sum(axis=1)
                   .fillna(1) > 0)]                      # drop consecutive duplicates
    out.to_csv(TAB / "stairway_ne_set1.csv", index=False)
    print(f"[31a] Stairway cleaned: {len(out)} steps "
          f"(year {out['year'].min():.1f}-{out['year'].max():.0f}); "
          f"recent Ne_median = {out['ne_median'].iloc[0]:.0f}")


# ---------- (b) within-largest-cluster D and Ne ----------
def largest_cluster_members():
    q = pd.read_csv(TAB / "admixture_set1_Q.csv")
    qcols = [c for c in q.columns if c.startswith("Q")]
    q["cl"] = "C" + q[qcols].to_numpy().argmax(axis=1).astype(str)
    top = q["cl"].value_counts().idxmax()
    return top, q.loc[q["cl"] == top, "sample"].astype(str).tolist()


def within_cluster_tajima(members):
    cs = allel.read_vcf(str(INTERIM / "Set1_nomaf.vcf"),
                        fields=["samples", "variants/CHROM", "variants/POS", "calldata/GT"])
    samples = list(cs["samples"])
    idx = [samples.index(m) for m in members if m in samples]
    gt = allel.GenotypeArray(cs["calldata/GT"][:, idx, :])
    chrom, pos = cs["variants/CHROM"], cs["variants/POS"]
    ac = gt.count_alleles()
    seg = ac.is_segregating()
    rows = []
    for ch in pd.unique(chrom):
        m = (chrom == ch) & seg
        if m.sum() < 50:
            continue
        pos_c, ac_c = pos[m], ac[m]
        order = np.argsort(pos_c, kind="stable")
        pos_c, ac_c = pos_c[order], ac_c[order]
        rows.append({"chrom": ch, "n_snps": int(m.sum()),
                     "tajima_d": allel.tajima_d(ac_c, pos=pos_c)})
    cw = pd.DataFrame(rows)
    w = cw["n_snps"] / cw["n_snps"].sum()
    return float((cw["tajima_d"] * w).sum()), int(cw["n_snps"].sum())


def within_cluster_ne(cluster, members, F_within):
    keep = INTERIM / "keep_largest_cluster.txt"
    keep.write_text("\n".join(f"{m} {m}" for m in members) + "\n")
    out = INTERIM / "ldwc"
    run_plink("--bfile", INTERIM / "Set1_qc", "--keep", keep, "--thin", "0.1",
              "--seed", "42", "--r2", "--ld-window-kb", "1000", "--ld-window", "99999",
              "--ld-window-r2", "0", "--out", out)
    ld = pd.read_csv(f"{out}.ld", sep=r"\s+")
    d = (ld["BP_B"] - ld["BP_A"]).abs()
    bins = (d // 25_000).astype(int)
    g = ld.groupby(bins)["R2"].mean().reset_index()
    g["bin_mid_bp"] = g["index"] * 25_000 + 12_500
    g = g[g["bin_mid_bp"] >= MIN_BIN_BP]
    n = len(members)
    rows = []
    for cm in CM_GRID:
        c = g["bin_mid_bp"].to_numpy() * cm * 1e-8
        r2adj = g["R2"].to_numpy() - 1.0 / n
        ok = r2adj > 0.01
        ne_pan = (1.0 / (4 * c[ok])) * (1.0 / r2adj[ok] - ALPHA)
        c_star = c[ok] * (1 - F_within)
        ne_adj = (1.0 / (4 * c_star)) * (1.0 / r2adj[ok] - ALPHA)
        order = np.argsort(-g["bin_mid_bp"].to_numpy()[ok])[:3]
        rec_pan = float(len(order) / (1.0 / ne_pan[order]).sum())
        rec_adj = float(len(order) / (1.0 / ne_adj[order]).sum())
        rows.append({"cm_per_mb": cm, "recent_ne_panmictic": rec_pan,
                     "recent_ne_selfing_adj": rec_adj})
    return pd.DataFrame(rows)


# ---------- (c) Ne interval across F choices ----------
def ne_interval():
    hwe = pd.read_csv(TAB / "hwe_summary.csv").set_index("panel")
    roh = pd.read_csv(TAB / "roh_summary.csv").set_index("panel")
    rows = []
    for panel in ["Set1", "Set2"]:
        ld = pd.read_csv(TAB / f"ld_decay_{panel.lower()}.csv")
        ld = ld[(ld["panel"] == panel) & (ld["bin_mid_bp"] >= MIN_BIN_BP)]
        n = {"Set1": 150, "Set2": 147}[panel]
        for F_name, F in [("PLINK_F", hwe.loc[panel, "mean_F"]),
                          ("F_ROH", roh.loc[panel, "mean_f_roh"])]:
            for cm in CM_GRID:
                c = ld["bin_mid_bp"].to_numpy() * cm * 1e-8
                r2adj = ld["mean_r2"].to_numpy() - 1.0 / n
                ok = r2adj > 0.01
                ne_pan = (1.0 / (4 * c[ok])) * (1.0 / r2adj[ok] - ALPHA)
                cs = c[ok] * (1 - F)
                ne_adj = (1.0 / (4 * cs)) * (1.0 / r2adj[ok] - ALPHA)
                order = np.argsort(-ld["bin_mid_bp"].to_numpy()[ok])[:3]
                rec_adj = float(len(order) / (1.0 / ne_adj[order]).sum())
                rec_pan = float(len(order) / (1.0 / ne_pan[order]).sum())
                t_epoch = float(np.mean(1.0 / (2 * cs[order])))
                rows.append({"panel": panel, "F_source": F_name, "F": float(F),
                             "cm_per_mb": cm, "recent_ne_panmictic": rec_pan,
                             "recent_ne_selfing_adj": rec_adj,
                             "epoch_generations_selfing_adj": t_epoch})
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "ne_interval_summary.csv", index=False)
    for panel in ["Set1", "Set2"]:
        sub = df[df["panel"] == panel]
        print(f"[31c] {panel}: Ne(selfing-adj) interval "
              f"{sub['recent_ne_selfing_adj'].min():.0f}-"
              f"{sub['recent_ne_selfing_adj'].max():.0f} "
              f"(panmictic {sub['recent_ne_panmictic'].min():.0f}-"
              f"{sub['recent_ne_panmictic'].max():.0f}); "
              f"epoch ~{sub['epoch_generations_selfing_adj'].mean():.0f} gen")
    return df


def main():
    clean_stairway()
    cluster, members = largest_cluster_members()
    fis = pd.read_csv(TAB / "fis_per_accession_set1.csv")
    F_within = float(fis[fis["sample"].isin(members)]["F"].mean())
    d_within, n_snps = within_cluster_tajima(members)
    print(f"[31b] largest cluster {cluster} (n={len(members)}, F={F_within:.3f}): "
          f"within-cluster Tajima's D = {d_within:+.2f} ({n_snps} segregating SNPs) "
          f"vs pooled +1.14")
    ne_wc = within_cluster_ne(cluster, members, F_within)
    ne_wc.insert(0, "cluster", cluster)
    ne_wc.insert(1, "n", len(members))
    base = pd.DataFrame([{"cluster": cluster, "n": len(members), "F_within": F_within,
                          "tajima_d_within": d_within, "n_snps_within": n_snps}])
    base.merge(ne_wc, on=["cluster", "n"]).to_csv(
        TAB / "within_cluster_history_set1.csv", index=False)
    ne_interval()
    print("[31] done")


if __name__ == "__main__":
    main()
