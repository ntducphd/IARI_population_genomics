#!/usr/bin/env python
"""28_figures_flagship.py — figures for the flagship-upgrade stages (13-26).

Produces (PNG + PDF via figstyle.save; panels skip gracefully when a stage has not run):
  main/Fig06_genome_history   (a) ROH length-class spectrum per panel
                              (b) F_ROH by admixture cluster (Set 1)
                              (c) LD-based Ne trajectory (panmictic vs selfing-adjusted)
                              (d) folded SFS, Set 1 (no-MAF, haploidised)
  main/Fig07_selection        (a,b) pcadapt Manhattan Set 1 / Set 2 (BH q<0.05 highlighted)
                              (c,d) iHS and XP-EHH Manhattans (Set 1) when stage 17 has run
  supp/SuppFig03_stairway     Stairway Plot 2 Ne trajectory (median + 95% CI) when stage 19 done
  supp/SuppFig04_f3           Patterson f3 Z-scores by target group (Z = -3 admixture line)
  supp/SuppFig05_dapc         DAPC vs admixture agreement (confusion heatmaps + ARI)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import TAB, FIG_MAIN, FIG_SUPP
import figstyle as fs

fs.set_style()
PANELS = ["Set1", "Set2"]
CHROM_SHADES = ["#4c4c4c", "#9a9a9a"]


def have(*names):
    return all((TAB / n).exists() for n in names)


def manhattan(ax, df, chrom_col, pos_col, y, highlight=None, ylabel=""):
    df = df.dropna(subset=[y]).copy()
    df[chrom_col] = df[chrom_col].astype(str).str.replace("chr", "", case=False).astype(int)
    df = df.sort_values([chrom_col, pos_col])
    offset, ticks = 0, []
    for i, (ch, sub) in enumerate(df.groupby(chrom_col)):
        x = sub[pos_col].to_numpy() + offset
        ax.scatter(x, sub[y], s=1.2, c=CHROM_SHADES[i % 2], rasterized=True, linewidths=0)
        if highlight is not None and highlight in sub.columns:
            hl = sub[sub[highlight]]
            ax.scatter(hl[pos_col].to_numpy() + offset, hl[y], s=2.5, c="#c0392b",
                       rasterized=True, linewidths=0)
        ticks.append(offset + sub[pos_col].max() / 2)
        offset += sub[pos_col].max() + 5e6
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(c) for c in sorted(df[chrom_col].unique())])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(ylabel)


def fig08():
    if not have("roh_length_classes_set1.csv", "roh_indiv_set1.csv", "ne_trajectory_set1.csv",
                "sfs_folded_set1.csv"):
        print("[28] Fig08 skipped (stages 13/14/19 outputs missing)")
        return
    fig, axes = plt.subplots(2, 2, figsize=(fs.DOUBLE, fs.DOUBLE * 0.62))
    (a, b), (c, d) = axes

    # (a) ROH length classes
    width = 0.38
    for i, p in enumerate(PANELS):
        lc = pd.read_csv(TAB / f"roh_length_classes_{p.lower()}.csv")
        x = np.arange(len(lc))
        a.bar(x + (i - 0.5) * width, 100 * lc["share_of_roh"], width,
              color=fs.PANEL_COL[p], label=fs.PANEL_LABEL[p])
        a.set_xticks(x)
        a.set_xticklabels(lc["length_class"], rotation=30, ha="right")
    a.set_ylabel("Share of total ROH length (%)")
    a.legend()
    fs.panel_letter(a, "a")

    # (b) F_ROH by cluster, Set 1
    indiv = pd.read_csv(TAB / "roh_indiv_set1.csv").dropna(subset=["cluster"])
    order = sorted(indiv["cluster"].unique())
    data = [indiv.loc[indiv["cluster"] == cl, "f_roh"] for cl in order]
    bp = b.boxplot(data, tick_labels=order, showfliers=False, widths=0.55,
                   patch_artist=True, medianprops={"color": "black"})
    for patch in bp["boxes"]:
        patch.set_facecolor(fs.PANEL_COL["Set1"])
        patch.set_alpha(0.45)
    for i, vals in enumerate(data):
        b.scatter(np.full(len(vals), i + 1) + np.random.default_rng(1).uniform(
            -0.12, 0.12, len(vals)), vals, s=3, c="#333333", alpha=0.6, linewidths=0)
    b.set_xlabel("Admixture cluster (Set 1)")
    b.set_ylabel("F_ROH")
    fs.panel_letter(b, "b")

    # (c) Ne trajectory
    for p in PANELS:
        tr = pd.read_csv(TAB / f"ne_trajectory_{p.lower()}.csv").dropna(
            subset=["ne_panmictic"])
        c.plot(tr["t_gen_panmictic"], tr["ne_panmictic"], color=fs.PANEL_COL[p],
               label=f"{fs.PANEL_LABEL[p]} (panmictic)")
        c.plot(tr["t_gen_selfing_adj"], tr["ne_selfing_adj"], color=fs.PANEL_COL[p],
               linestyle="--", label=f"{fs.PANEL_LABEL[p]} (selfing-adj.)")
    c.set_xscale("log"); c.set_yscale("log")
    c.set_xlabel("Generations ago (~ 1/(2c))")
    c.set_ylabel("Ne (LD-based)")
    c.legend(fontsize=5.6)
    # the deep-past rise is a Sved-binning mechanical artefact the caption
    # tells readers not to interpret — say so ON the figure, in the region it applies to
    c.axvspan(200, c.get_xlim()[1], color="0.92", zorder=0)
    c.annotate("mechanical rise\n(not interpreted)", xy=(0.97, 0.06),
               xycoords="axes fraction", ha="right", va="bottom",
               fontsize=fs.FS_ANNOT, color="0.35")
    fs.panel_letter(c, "c")

    # (d) folded SFS
    sfs = pd.read_csv(TAB / "sfs_folded_set1.csv")
    d.bar(sfs["minor_allele_count"], sfs["n_snps"], width=1.0,
          color=fs.PANEL_COL["Set1"])
    d.set_xlabel("Minor-allele count (n = 150 haploidised)")
    d.set_ylabel("SNPs")
    d.set_yscale("log")
    fs.panel_letter(d, "d")

    fig.tight_layout()
    fs.save(fig, str(FIG_MAIN / "Fig06_genome_history.png"))
    plt.close(fig)
    print("[28] Fig06_genome_history written")


def fig09():
    if not have("pcadapt_outliers_set1.csv", "pcadapt_outliers_set2.csv"):
        print("[28] Fig09 skipped (stage 18 missing)")
        return
    has17 = have("ihs_set1.csv", "xpehh_set1.csv")
    n_rows = 2 if has17 else 1
    fig, axes = plt.subplots(n_rows, 2, figsize=(fs.DOUBLE, fs.DOUBLE * 0.33 * n_rows),
                             squeeze=False)

    for i, p in enumerate(PANELS):
        pc = pd.read_csv(TAB / f"pcadapt_outliers_{p.lower()}.csv")
        pc["logp"] = -np.log10(pc["pvalue"].clip(lower=1e-300))
        manhattan(axes[0][i], pc, "chrom", "pos", "logp", highlight="outlier",
                  ylabel="-log10 P (pcadapt)")
        axes[0][i].set_title(fs.PANEL_LABEL[p])
        if p == "Set2":
            # the visually dramatic band is the numerical floor of the
            # test, not evidence strength — annotate it where the eye lands
            top = pc["logp"].max()
            axes[0][i].annotate("numerical P floor (saturated) — see text",
                               xy=(0.5, top), xycoords=("axes fraction", "data"),
                               ha="center", va="bottom", fontsize=fs.FS_ANNOT,
                               color="0.35")
        fs.panel_letter(axes[0][i], "ab"[i])

    if has17:
        ihs = pd.read_csv(TAB / "ihs_set1.csv")
        ihs["abs_ihs"] = ihs["ihs"].abs()
        manhattan(axes[1][0], ihs, "chrom", "pos", "abs_ihs", highlight="outlier",
                  ylabel="|iHS| (Set 1)")
        axes[1][0].axhline(3, lw=0.6, ls=":", color="#c0392b")
        fs.panel_letter(axes[1][0], "c")
        xp = pd.read_csv(TAB / "xpehh_set1.csv")
        xp["abs_xp"] = xp["xpehh"].abs()
        manhattan(axes[1][1], xp, "chrom", "pos", "abs_xp", highlight="outlier",
                  ylabel="|XP-EHH| (Set 1)")
        axes[1][1].axhline(3, lw=0.6, ls=":", color="#c0392b")
        fs.panel_letter(axes[1][1], "d")

    fig.tight_layout()
    fs.save(fig, str(FIG_MAIN / "Fig07_selection.png"))
    plt.close(fig)
    print(f"[28] Fig07_selection written ({'4' if has17 else '2'} panels)")


def supp_stairway():
    if not have("stairway_ne_set1.csv"):
        print("[28] SuppFig04 skipped (stage 19 pending)")
        return
    st = pd.read_csv(TAB / "stairway_ne_set1.csv")
    # belt-and-braces: exclude numeric-underflow rows even if the table was not pre-cleaned
    st = st[np.isfinite(st["year"]) & (st["year"] >= 0.1)]
    fig, ax = plt.subplots(figsize=(fs.ONE_HALF, fs.ONE_HALF * 0.7))
    ax.plot(st["year"], st["ne_median"], color=fs.PANEL_COL["Set1"], label="Ne median")
    if {"ne_lo95", "ne_hi95"}.issubset(st.columns):
        ax.fill_between(st["year"], st["ne_lo95"], st["ne_hi95"],
                        color=fs.PANEL_COL["Set1"], alpha=0.2, linewidth=0,
                        label="95% CI (exploratory)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Years ago (mu = 7e-9, 1 gen/yr — scaling constants)")
    ax.set_ylabel("Ne (Stairway Plot 2)")
    ax.legend()
    fig.tight_layout()
    fs.save(fig, str(FIG_SUPP / "SuppFig03_stairway.png"))
    plt.close(fig)
    print("[28] SuppFig03_stairway written")


def supp_f3():
    if not have("f3_set1.csv"):
        return
    f3 = pd.read_csv(TAB / "f3_set1.csv")
    targets = sorted(f3["target"].unique())
    fig, ax = plt.subplots(figsize=(fs.ONE_HALF, fs.ONE_HALF * 0.7))
    rng = np.random.default_rng(2)
    for i, t in enumerate(targets):
        z = f3.loc[f3["target"] == t, "z"]
        ax.scatter(np.full(len(z), i) + rng.uniform(-0.15, 0.15, len(z)), z,
                   s=14, c=np.where(z <= -3, "#c0392b", "#666666"), linewidths=0)
    ax.axhline(-3, lw=0.7, ls="--", color="#c0392b")
    ax.annotate("Z = -3 (admixed)", xy=(0.99, -3), xycoords=("axes fraction", "data"),
                ha="right", va="bottom", fontsize=fs.FS_ANNOT, color="#c0392b")
    ax.set_xticks(range(len(targets)))
    ax.set_xticklabels(targets)
    ax.set_xlabel("Target group (3K-RGP labels, Set 1)")
    ax.set_ylabel("Patterson f3 Z-score")
    fig.tight_layout()
    fs.save(fig, str(FIG_SUPP / "SuppFig04_f3.png"))
    plt.close(fig)
    print("[28] SuppFig04_f3 written")


def supp_dapc():
    if not have("dapc_assign_set1.csv", "dapc_assign_set2.csv"):
        return
    fig, axes = plt.subplots(1, 2, figsize=(fs.DOUBLE, fs.DOUBLE * 0.42))
    for ax, p in zip(axes, PANELS):
        dapc = pd.read_csv(TAB / f"dapc_assign_{p.lower()}.csv")
        q = pd.read_csv(TAB / f"admixture_{p.lower()}_Q.csv")
        qcols = [c for c in q.columns if c.startswith("Q")]
        q["admix"] = "C" + q[qcols].to_numpy().argmax(axis=1).astype(str)
        m = dapc.merge(q[["sample", "admix"]], on="sample")
        ari = adjusted_rand_score(m["admix"], m["dapc_cluster"])
        ct = pd.crosstab(m["admix"], m["dapc_cluster"])
        im = ax.imshow(ct.to_numpy(), cmap="Blues", aspect="auto")
        ax.set_xticks(range(len(ct.columns))); ax.set_xticklabels(ct.columns, fontsize=5.5)
        ax.set_yticks(range(len(ct.index))); ax.set_yticklabels(ct.index, fontsize=5.5)
        ax.set_xlabel("DAPC cluster"); ax.set_ylabel("Admixture cluster")
        ax.set_title(f"{fs.PANEL_LABEL[p]} — ARI = {ari:.2f}")
        for (i, j), v in np.ndenumerate(ct.to_numpy()):
            if v:
                ax.annotate(str(v), (j, i), ha="center", va="center", fontsize=5,
                            color="white" if v > ct.to_numpy().max() / 2 else "#222222")
    fig.tight_layout()
    fs.save(fig, str(FIG_SUPP / "SuppFig05_dapc.png"))
    plt.close(fig)
    print("[28] SuppFig05_dapc written")


if __name__ == "__main__":
    fig08()
    fig09()
    supp_stairway()
    supp_f3()
    supp_dapc()
    print("[28] done")

