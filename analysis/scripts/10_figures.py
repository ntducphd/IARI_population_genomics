#!/usr/bin/env python
"""10_figures.py — the stage-10 figures in the final high-density figure scheme.
PNG for review + vector PDF for the journal (figstyle.save).

Fig 3   NJ tree + GRM kinship heatmap, both panels (tips coloured by admixture cluster)
Supp S1 UMAP embedding + structure-consensus (ARI/NMI) heatmap
Supp S2 3K-RGP global anchor PCA (Set1 only; Set2 gap noted)

The remaining floats in the final scheme (structure, differentiation/diversity/LD, concordance,
core-collection application) render as composites in 35_figures_composite.py; the design schematic
is built by 29_fig01_design.py.
"""
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as fs
from paths import TAB, INTERIM, FIG_MAIN, FIG_SUPP

fs.set_style()
PANELS = ["Set1", "Set2"]


def cluster_palette(n):
    # shared Okabe-Ito colour-blind-safe cluster palette (figstyle.CLUSTER_COL), indexed C0..;
    # replaced the matplotlib tab10 default 2026-08-08 (figure-standard pass)
    return [fs.cluster_color(i) for i in range(n)]


# 3K-RGP anchor figure needs MORE colours than the 9-colour cluster palette (11 labelled
# subpopulation groups) — cycling the palette silently duplicated aus/trop and admix/temp
# (QA catch 2026-08-08). Two extra distinguishable colours close the gap.
ANCHOR_EXTRA = ["#8c510a", "#01665e"]


def anchor_palette(n):
    base = [fs.cluster_color(i) for i in range(min(n, 9))]
    return (base + ANCHOR_EXTRA)[:n]


def method_label(name):
    # structure-consensus method codes -> reader-facing labels (QA pass 2026-08-08)
    if name == "admixture_argmaxQ":
        return "Admixture\n(dominant Q)"
    m = re.match(r"(PCA|UMAP)_kmeans_k(\d+)", name)
    return f"{m.group(1)} k-means\n(k = {m.group(2)})" if m else name


# =============================================================================
# Fig 3 — NJ tree + kinship heatmap (final-scheme name: Fig03_tree_kinship)
# =============================================================================
def fig3_tree_kinship():
    fig = plt.figure(figsize=(fs.DOUBLE, fs.DOUBLE * 0.95))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.3, 1], hspace=0.35, wspace=0.30)
    for i, panel in enumerate(PANELS):
        ax = fig.add_subplot(gs[0, i])
        tree = Phylo.read(str(INTERIM / f"nj_{panel.lower()}.nwk"), "newick")
        # tip labels suppressed via label_func (names must be KEPT for cluster colouring below)
        Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False,
                    label_func=lambda c: "")
        # colour tip markers by admixture cluster (same palette as Fig 2 — the teaching link)
        q = pd.read_csv(TAB / f"admixture_{panel.lower()}_Q.csv")
        qcols = [c for c in q.columns if c.startswith("Q")]
        clmap = dict(zip(q["sample"].astype(str),
                         q[qcols].to_numpy().argmax(axis=1)))
        depths = tree.depths()
        for yi, term in enumerate(tree.get_terminals(), start=1):
            cl = clmap.get(str(term.name))
            if cl is not None:
                ax.scatter(depths[term], yi, s=2.6, color=fs.cluster_color(int(cl)),
                           zorder=4, linewidths=0)
        ax.set_yticks([])                       # tip indices are meaningless — hide them
        ax.set_title(f"{fs.PANEL_LABEL[panel]} NJ tree ({tree.count_terminals()} tips)",
                     loc="left", fontsize=fs.FS_TITLE - 0.5)
        ax.set_xlabel("IBS distance"); ax.set_ylabel("")
        fs.panel_letter(ax, "ab"[i])

        ax2 = fig.add_subplot(gs[1, i])
        kin = pd.read_csv(TAB / f"kinship_{panel.lower()}.csv")
        ids = sorted(set(kin["sample_a"]))
        k = kin.pivot(index="sample_a", columns="sample_b", values="kinship").reindex(index=ids, columns=ids)
        im = ax2.imshow(k.to_numpy(), cmap="viridis", aspect="auto")
        ax2.set_xticks([]); ax2.set_yticks([])
        ax2.set_title(f"{fs.PANEL_LABEL[panel]} GRM kinship ({len(ids)}x{len(ids)})", loc="left",
                       fontsize=fs.FS_TITLE - 0.5)
        cbar = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=fs.FS_TICK - 1)
        fs.panel_letter(ax2, "cd"[i])
    fs.save(fig, str(FIG_MAIN / "Fig03_tree_kinship.png"))
    plt.close(fig)
    print("  Fig 3 (tree + kinship) done")


# =============================================================================
# Supp Fig S1 — UMAP + structure consensus (ARI/NMI) heatmap
# =============================================================================
def suppfig1_umap_consensus():
    fig, axes = plt.subplots(2, 2, figsize=(fs.DOUBLE, fs.DOUBLE * 0.95))
    fig.subplots_adjust(hspace=0.55, wspace=0.3)
    for i, panel in enumerate(PANELS):
        ax = axes[0, i]
        um = pd.read_csv(TAB / f"umap_{panel.lower()}.csv")
        q = pd.read_csv(TAB / f"admixture_{panel.lower()}_Q.csv")
        qcols = [c for c in q.columns if c.startswith("Q")]
        q["cluster"] = q[qcols].to_numpy().argmax(axis=1)
        m = um.merge(q[["sample", "cluster"]], on="sample")
        cols = cluster_palette(m["cluster"].nunique())
        for k, cl in enumerate(sorted(m["cluster"].unique())):
            sub = m[m["cluster"] == cl]
            ax.scatter(sub["UMAP1"], sub["UMAP2"], s=12, color=cols[k], edgecolor="0.2", linewidth=0.2)
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
        ax.set_title(f"{fs.PANEL_LABEL[panel]} UMAP (coloured by admixture cluster)", loc="left",
                     fontsize=fs.FS_TITLE - 1)
        fs.panel_letter(ax, "ab"[i])

        ax2 = axes[1, i]
        sc = pd.read_csv(TAB / f"structure_consensus_{panel.lower()}.csv")
        methods = sorted(set(sc["method_a"]) | set(sc["method_b"]))
        M = np.full((len(methods), len(methods)), np.nan)
        np.fill_diagonal(M, 1.0)
        for _, r in sc.iterrows():
            a, b = methods.index(r.method_a), methods.index(r.method_b)
            M[a, b] = M[b, a] = r.ARI
        im = ax2.imshow(np.ma.masked_invalid(M), cmap="Greens", vmin=0, vmax=1)
        im.cmap.set_bad("0.92")            # not-computed pairs: grey, annotated below
        ax2.set_xticks(range(len(methods)))
        ax2.set_xticklabels([method_label(m) for m in methods], rotation=45, fontsize=5,
                            ha="right")
        ax2.set_yticks(range(len(methods)))
        ax2.set_yticklabels([method_label(m) for m in methods], fontsize=5)
        for a in range(len(methods)):
            for b in range(len(methods)):
                if np.isfinite(M[a, b]):
                    ax2.text(b, a, f"{M[a,b]:.2f}", ha="center", va="center", fontsize=5,
                             color="white" if M[a, b] > 0.5 else "0.1")
                else:
                    ax2.text(b, a, "–", ha="center", va="center", fontsize=6, color="0.45")
        ax2.set_title(f"{fs.PANEL_LABEL[panel]} cross-method ARI", loc="left", fontsize=fs.FS_TITLE - 1)
        fs.panel_letter(ax2, "cd"[i])
    fs.save(fig, str(FIG_SUPP / "SuppFig01_umap_consensus.png"))
    plt.close(fig)
    print("  Supp Fig S1 done")


# =============================================================================
# Supp Fig S2 — 3K-RGP global anchor
# =============================================================================
def suppfig2_anchor():
    fig, ax = plt.subplots(figsize=(fs.ONE_HALF, fs.ONE_HALF * 0.85))
    df = pd.read_csv(TAB / "anchor_pca_global.csv")
    var = pd.read_csv(TAB / "anchor_pca_variance.csv").set_index("PC")
    ref = df[df["group"] == "3K-RGP reference (unlabelled)"]
    ax.scatter(ref["PC1"], ref["PC2"], s=3, color="0.82", label=f"3K-RGP reference (n={len(ref)})", zorder=1)
    labelled = df[df["group"] != "3K-RGP reference (unlabelled)"]
    groups = sorted(labelled["group"].unique())
    cols = anchor_palette(len(groups))
    for k, g in enumerate(groups):
        sub = labelled[labelled["group"] == g]
        ax.scatter(sub["PC1"], sub["PC2"], s=16, color=cols[k], edgecolor="0.15", linewidth=0.3,
                   label=f"{g} (n={len(sub)})", zorder=3)
    ax.set_xlabel(f"PC1 ({var.loc['PC1','pct_variance']:.1f}%)")
    ax.set_ylabel(f"PC2 ({var.loc['PC2','pct_variance']:.1f}%)")
    # single-panel figure: an in-image title would be the figure title — banned (QA 2026-08-08)
    ax.legend(fontsize=4.4, ncol=2, loc="best")
    fs.save(fig, str(FIG_SUPP / "SuppFig02_anchor_3krgp.png"))
    plt.close(fig)
    print("  Supp Fig S2 done")


if __name__ == "__main__":
    print("=== Building stage-10 figures (final 8-figure scheme) ===")
    fig3_tree_kinship()
    suppfig1_umap_consensus()
    suppfig2_anchor()
    print(f"\n-> {FIG_MAIN} (Fig03_tree_kinship)")
    print(f"-> {FIG_SUPP} (SuppFig01, SuppFig02)")
