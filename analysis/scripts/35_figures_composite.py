#!/usr/bin/env python
"""35_figures_composite.py — [stage 35] high-density figure consolidation (design rule:
4-8 panels per main figure). FINAL 10-figure scheme (2026-08-08, author-confirmed — the two
extra mains are genuinely NEW figures, not supp promotions):

  Fig01 design (29_fig01_design.py)          Fig06 genome history (28)
  Fig02 structure & ancestry (THIS, 6)       Fig07 selection (28)
  Fig03 trees & kinship (10_figures.py)      Fig08 application (THIS, 3)
  Fig04 differentiation/diversity/LD (THIS)  Fig09 mechanism made visible (THIS, 6) <- stage 36
  Fig05 concordance+mechanism (THIS, 6)      Fig10 two platforms, one biology (THIS, 3) <- 36

Shared cluster palette: figstyle.CLUSTER_COL (C0.. labels). No in-image figure titles.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import TAB, FIG_MAIN
import figstyle as fs

fs.set_style()
PANELS = ["Set1", "Set2"]


def have(*names):
    missing = [n for n in names if not (TAB / n).exists()]
    if missing:
        print(f"  (missing {missing})")
    return not missing


# ============================================================================
# Fig 2 — structure & ancestry (6 panels)
# ============================================================================
def fig02_structure():
    fig = plt.figure(figsize=(fs.DOUBLE, fs.DOUBLE * 0.92))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 0.9, 0.9], hspace=0.55, wspace=0.28)
    var = pd.read_csv(TAB / "pca_variance.csv").set_index("panel")
    bestk = pd.read_csv(TAB / "admixture_bestK.csv").set_index("panel")

    for i, panel in enumerate(PANELS):
        # (a,b) PCA coloured by cluster
        ax = fig.add_subplot(gs[0, i])
        pca = pd.read_csv(TAB / f"pca_{panel.lower()}.csv")
        q = pd.read_csv(TAB / f"admixture_{panel.lower()}_Q.csv")
        qcols = [c for c in q.columns if c.startswith("Q")]
        q["cluster"] = q[qcols].to_numpy().argmax(axis=1)
        m = pca.merge(q[["sample", "cluster"]], on="sample", how="inner")
        for cl in sorted(m["cluster"].unique()):
            sub = m[m["cluster"] == cl]
            ax.scatter(sub["PC1"], sub["PC2"], s=13, color=fs.cluster_color(int(cl)),
                       edgecolor="0.2", linewidth=0.25, label=fs.cluster_label(int(cl)),
                       zorder=3)
        ax.set_xlabel(f"PC1 ({var.loc[panel, 'PC1']:.1f}%)")
        ax.set_ylabel(f"PC2 ({var.loc[panel, 'PC2']:.1f}%)")
        ax.set_title(f"{fs.PANEL_LABEL[panel]}, n = {len(m)}", loc="left",
                     fontsize=fs.FS_TITLE - 0.5)
        # deterministic legend spots verified point-free (QA sweep): a upper-left, b upper-centre
        ax.legend(fontsize=4.4, ncol=2, loc="upper left" if i == 0 else "upper center",
                  handletextpad=0.3, borderaxespad=0.2)
        fs.panel_letter(ax, "ab"[i])

        # (c,d) cross-entropy K-selection
        ax2 = fig.add_subplot(gs[1, i])
        ce = pd.read_csv(TAB / f"admixture_{panel.lower()}_cv.csv")
        mean_ce = ce.groupby("K")["cross_entropy"].agg(["mean", "std"]).reset_index()
        ax2.errorbar(mean_ce["K"], mean_ce["mean"], yerr=mean_ce["std"],
                     color=fs.PANEL_COL[panel], marker="o", ms=3, lw=1.0, capsize=2)
        elbow = int(bestk.loc[panel, "bestK_elbow"])
        gmin = int(bestk.loc[panel, "bestK_global_min"])
        ax2.axvline(elbow, color="0.3", ls="--", lw=0.8)
        ax2.text(elbow - 0.15, ax2.get_ylim()[1], f"elbow K={elbow} ", fontsize=5.6,
                 va="top", ha="right", color="0.2")
        ax2.scatter([gmin], [mean_ce.loc[mean_ce.K == gmin, "mean"].iloc[0]], marker="*",
                    s=40, color="firebrick", zorder=4, label=f"global-min K={gmin}")
        ax2.set_xlabel("K (ancestral populations)")
        ax2.set_ylabel("Cross-entropy")
        ax2.legend(fontsize=5.2, loc="lower left")
        fs.panel_letter(ax2, "cd"[i])

        # (e,f) Q-matrix stacked bars
        ax3 = fig.add_subplot(gs[2, i])
        dom = q[qcols].to_numpy().argmax(axis=1)
        order = np.lexsort((-q[qcols].to_numpy().max(axis=1), dom))
        qs = q.iloc[order].reset_index(drop=True)
        bottom = np.zeros(len(qs))
        for k, qc in enumerate(qcols):
            ax3.bar(np.arange(len(qs)), qs[qc], bottom=bottom, width=1.0,
                    color=fs.cluster_color(k), linewidth=0)
            bottom += qs[qc].to_numpy()
        ax3.set_xlim(-0.5, len(qs) - 0.5)
        ax3.set_ylim(0, 1)
        ax3.set_xticks([])
        ax3.set_xlabel("Accessions (sorted by dominant ancestry)")
        ax3.set_ylabel("Ancestry proportion")
        ax3.set_title(f"K = {len(qcols)}", loc="left", fontsize=fs.FS_TITLE - 0.5)
        fs.panel_letter(ax3, "ef"[i])

    fs.save(fig, str(FIG_MAIN / "Fig02_structure.png"))
    plt.close(fig)
    print("[35] Fig02_structure written (6 panels)")


# ============================================================================
# Fig 4 — differentiation, diversity, LD (4 panels)
# ============================================================================
def fig04_diff_diversity_ld():
    fig, axes = plt.subplots(2, 2, figsize=(fs.DOUBLE, fs.DOUBLE * 0.80),
                             gridspec_kw={"hspace": 0.45, "wspace": 0.35})
    (a, b), (c, d) = axes

    for ax, panel, letter in zip((a, b), PANELS, "ab"):
        pw = pd.read_csv(TAB / f"fst_pairwise_{panel.lower()}.csv")
        clusters = sorted(set(pw["cluster_a"]) | set(pw["cluster_b"]))
        M = np.zeros((len(clusters), len(clusters)))
        for _, r in pw.iterrows():
            ia, ib = clusters.index(r.cluster_a), clusters.index(r.cluster_b)
            M[ia, ib] = M[ib, ia] = r.fst
        im = ax.imshow(M, cmap="viridis", vmin=0)
        labels = [fs.cluster_label(int(cc)) for cc in clusters]
        ax.set_xticks(range(len(clusters))); ax.set_xticklabels(labels, fontsize=5.0)
        ax.set_yticks(range(len(clusters))); ax.set_yticklabels(labels, fontsize=5.0)
        vmax = M.max() if M.max() > 0 else 1.0
        for ia in range(len(clusters)):
            for ib in range(len(clusters)):
                if ia != ib:
                    ax.annotate(f"{M[ia, ib]:.2f}".lstrip("0"), (ib, ia), ha="center",
                                va="center", fontsize=3.4,
                                color="white" if M[ia, ib] < 0.55 * vmax else "black")
        glob = pd.read_csv(TAB / f"fst_global_{panel.lower()}.csv").iloc[0]
        ax.set_title(f"{fs.PANEL_LABEL[panel]} — Gst = {glob.global_fst:.3f} "
                     f"(P = {glob.perm_p:.3f})", loc="left", fontsize=fs.FS_TITLE - 1)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).ax.tick_params(
            labelsize=fs.FS_TICK - 1)
        fs.panel_letter(ax, letter)

    # (c) Ho/He/F
    div = pd.read_csv(TAB / "diversity_summary.csv")
    metrics = ["mean_Ho", "mean_He", "mean_F"]
    x = np.arange(len(metrics)); w = 0.35
    for i, panel in enumerate(PANELS):
        row = div[div.panel == panel].iloc[0]
        c.bar(x + (i - 0.5) * w, [row[mm] for mm in metrics], width=w,
              color=fs.PANEL_COL[panel], label=fs.PANEL_LABEL[panel])
    c.set_xticks(x); c.set_xticklabels(["Ho", "He", "F"])
    c.set_ylabel("Estimate")
    c.legend(fontsize=5.4)
    fs.panel_letter(c, "c")

    # (d) LD decay
    summ = pd.read_csv(TAB / "ld_decay_summary.csv").set_index("panel")
    for panel in PANELS:
        dd = pd.read_csv(TAB / f"ld_decay_{panel.lower()}.csv")
        d.plot(dd["bin_mid_bp"] / 1000, dd["mean_r2"], color=fs.PANEL_COL[panel],
               marker="o", ms=2.2, lw=1.0,
               label=f"{fs.PANEL_LABEL[panel]} (half-decay "
                     f"{summ.loc[panel, 'half_decay_bp'] / 1000:.0f} kb)")
        d.axvline(summ.loc[panel, "half_decay_bp"] / 1000, color=fs.PANEL_COL[panel],
                  ls=":", lw=0.8)
    d.axhline(0.2, color="#888888", ls="--", lw=0.7)
    d.annotate(r"$r^2$ = 0.2", xy=(0.985, 0.2), xycoords=("axes fraction", "data"),
               ha="right", va="bottom", fontsize=fs.FS_ANNOT, color="#555555")
    d.set_xlabel("Physical distance (kb)")
    d.set_ylabel(r"Mean $r^2$")
    d.legend(fontsize=5.2)
    fs.panel_letter(d, "d")

    fig.tight_layout()
    fs.save(fig, str(FIG_MAIN / "Fig04_diff_diversity_ld.png"))
    plt.close(fig)
    print("[35] Fig04_diff_diversity_ld written (4 panels)")


# ============================================================================
# Fig 5 — concordance + mechanism (6 panels; needs stage 30 complete)
# ============================================================================
def fig05_concordance_mechanism():
    need = ["concordance_mantel_set1.csv", "mantel_bootstrap_set1.csv",
            "mantel_bootstrap_set2.csv", "robustness_procrustes_set1.csv",
            "robustness_procrustes_set2.csv", "classifier_baselines_set1.csv",
            "classifier_baselines_set2.csv", "confound_mantel_set1.csv",
            "confound_mantel_set2.csv", "concordance_feature_attribution_set1.csv"]
    if not have(*need):
        print("[35] Fig05 skipped (stage 30 pending)")
        return
    fig, axes = plt.subplots(2, 3, figsize=(fs.DOUBLE, fs.DOUBLE * 0.66))
    (a, b, c), (d, e, f_) = axes

    # (a) Mantel r + bootstrap CI
    for i, p in enumerate(PANELS):
        m = pd.read_csv(TAB / f"concordance_mantel_{p.lower()}.csv")
        r = float(m.loc[m["comparison"] == "genomic~phenomic", "r"].iloc[0])
        bt = pd.read_csv(TAB / f"mantel_bootstrap_{p.lower()}.csv").iloc[0]
        a.errorbar(r, i, xerr=[[r - bt["ci95_lo"]], [bt["ci95_hi"] - r]],
                   fmt=fs.PANEL_MARKER[p], color=fs.PANEL_COL[p], ms=6, capsize=3, lw=1.2)
    a.axvline(0, lw=0.7, color="#888888", ls=":")
    a.set_yticks([0, 1]); a.set_yticklabels([fs.PANEL_LABEL[p] for p in PANELS])
    a.set_xlabel("Mantel r, genomic ~ phenomic\n(bootstrap 95% CI)")
    a.set_ylim(-0.7, 1.7)
    fs.panel_letter(a, "a")

    # (b) Procrustes band
    for p in PANELS:
        pr = pd.read_csv(TAB / f"robustness_procrustes_{p.lower()}.csv")
        b.plot(pr["n_pc"], pr["m2"], marker=fs.PANEL_MARKER[p], ms=4,
               color=fs.PANEL_COL[p], label=fs.PANEL_LABEL[p])
        prim = pr.loc[pr["n_pc"] == 4]
        b.scatter(prim["n_pc"], prim["m2"], s=70, facecolors="none",
                  edgecolors=fs.PANEL_COL[p], linewidths=1.4, zorder=5)
    b.set_xlabel("Principal components retained")
    b.set_ylabel("Procrustes M²")
    b.set_ylim(0.75, 1.0)
    b.legend(fontsize=5.4)
    fs.panel_letter(b, "b")

    # (c) classifier vs null + external labels
    xpos, w = np.arange(len(PANELS)), 0.35
    for i, p in enumerate(PANELS):
        cb = pd.read_csv(TAB / f"classifier_baselines_{p.lower()}.csv")
        row = cb[(cb["features"] == "full_204") &
                 (cb["target"] == "admixture_elbowK")].iloc[0]
        c.bar(i - w / 2, row["accuracy"], w, color=fs.PANEL_COL[p],
              label=None if i else "observed")
        c.bar(i + w / 2, row["null_mean"], w, color="#b9c2cc",
              label=None if i else "permutation null (mean)")
        c.annotate(f"P = {row['empirical_p']:.3f}", (i, row["accuracy"] + 0.015),
                   ha="center", fontsize=fs.FS_ANNOT)
    ext_p = TAB / "external_label_classifier_set1.csv"
    if ext_p.exists():
        e_row = pd.read_csv(ext_p).iloc[0]
        c.scatter([-w / 2], [e_row["accuracy"]], marker="D", s=26, color="#0f5f59",
                  zorder=5, label="external 3K-RGP labels")
    c.set_xticks(xpos); c.set_xticklabels([fs.PANEL_LABEL[p] for p in PANELS], fontsize=5.8)
    c.set_ylabel("Ancestry-classification accuracy")
    c.legend(fontsize=5.2, loc="upper right")
    fs.panel_letter(c, "c")

    # (d) family Mantel vs full
    fam_order = ["Size_morphology", "Colour", "NIR"]
    xpos = np.arange(len(fam_order))
    for i, p in enumerate(PANELS):
        fa = pd.read_csv(TAB / f"concordance_feature_attribution_{p.lower()}.csv")
        gcol = [cn for cn in fa.columns if "feat" in cn.lower() or "group" in cn.lower()][0]
        rcol = [cn for cn in fa.columns if cn.lower() in ("r", "mantel_r", "r_genomic")][0]
        vals = [float(fa.loc[fa[gcol] == fam, rcol].iloc[0]) for fam in fam_order]
        d.bar(xpos + (i - 0.5) * w, vals, w, color=fs.PANEL_COL[p],
              label=fs.PANEL_LABEL[p])
        m = pd.read_csv(TAB / f"concordance_mantel_{p.lower()}.csv")
        rfull = float(m.loc[m["comparison"] == "genomic~phenomic", "r"].iloc[0])
        d.axhline(rfull, lw=0.8, ls="--", color=fs.PANEL_COL[p], alpha=0.8)
    d.set_xticks(xpos)
    d.set_xticklabels(["Size/morph.", "Colour", "NIR"], fontsize=5.8)
    d.set_ylabel("Family Mantel r (dashed = full 204)")
    # opaque frame: the full-width dashed reference lines otherwise strike through the text
    d.legend(fontsize=5.2, frameon=True, facecolor="white", edgecolor="none", framealpha=1)
    fs.panel_letter(d, "d")

    # (e) stature-collapse dumbbell
    labels = ["full 204", "height only", "full | size\n(partial)", "non-size | size\n(partial)"]
    keymap = ["genomic~phenomic (primary)", "genomic~height_only",
              "genomic~phenomic | size_family", "genomic~nonsize | size_family"]
    ys = np.arange(len(labels))[::-1]
    for p in PANELS:
        cm = pd.read_csv(TAB / f"confound_mantel_{p.lower()}.csv").set_index("test")
        if not set(keymap).issubset(cm.index):
            continue
        vals = [float(cm.loc[k, "r"]) for k in keymap]
        ps = [float(cm.loc[k, "p"]) for k in keymap]
        e.plot(vals, ys, marker=fs.PANEL_MARKER[p], ms=5, lw=1.0,
               color=fs.PANEL_COL[p], label=fs.PANEL_LABEL[p])
        for v, y, pv in zip(vals, ys, ps):
            if pv >= 0.05:
                e.annotate("n.s.", (v, y + 0.16), ha="center", fontsize=5.4,
                           color=fs.PANEL_COL[p])
    e.axvline(0, lw=0.7, color="#888888", ls=":")
    e.set_yticks(ys); e.set_yticklabels(labels, fontsize=5.6)
    e.set_xlabel("Mantel r with genomic distance")
    e.legend(fontsize=5.2)
    fs.panel_letter(e, "e")

    # (f) classifier feature subsets
    subsets = ["full_204", "height_only", "nonsize_colour_nir"]
    slabels = ["full 204", "height only", "colour + NIR"]
    xpos = np.arange(len(subsets))
    for i, p in enumerate(PANELS):
        cb = pd.read_csv(TAB / f"classifier_baselines_{p.lower()}.csv")
        cb = cb[cb["target"] == "admixture_elbowK"].set_index("features")
        vals = [float(cb.loc[s, "accuracy"]) for s in subsets if s in cb.index]
        f_.bar(xpos[:len(vals)] + (i - 0.5) * w, vals, w, color=fs.PANEL_COL[p],
               label=fs.PANEL_LABEL[p])
        nm = cb.loc["full_204", "null_mean"]
        if pd.notna(nm):
            f_.axhline(float(nm), lw=0.8, ls=":", color=fs.PANEL_COL[p], alpha=0.8)
    f_.set_xticks(xpos); f_.set_xticklabels(slabels, fontsize=5.6)
    f_.set_ylabel("Accuracy (dotted = null)")
    # no legend: every candidate spot collides with a bar (QA sweep); panels b/d/e already
    # carry the shared Set1/Set2 legend for this figure
    fs.panel_letter(f_, "f")

    fig.tight_layout()
    fs.save(fig, str(FIG_MAIN / "Fig05_concordance_mechanism.png"))
    plt.close(fig)
    print("[35] Fig05_concordance_mechanism written (6 panels)")


# ============================================================================
# Fig 8 — application (3 panels)
# ============================================================================
def fig08_application():
    if not have("core_collection_curve_set1.csv", "richness_set1.csv"):
        print("[35] Fig08 skipped")
        return
    fig, axes = plt.subplots(1, 3, figsize=(fs.DOUBLE, fs.DOUBLE * 0.34))
    a, b, c = axes

    for p in PANELS:
        cv = pd.read_csv(TAB / f"core_collection_curve_{p.lower()}.csv")
        a.plot(cv["pct_of_panel"], cv["pct_diversity_retained"], marker=fs.PANEL_MARKER[p],
               ms=3, color=fs.PANEL_COL[p], label=fs.PANEL_LABEL[p])
    a.axvline(10, lw=0.8, ls="--", color="#888888")
    a.axhline(100, lw=0.7, ls=":", color="#888888")
    a.set_xlabel("Core size (% of panel)")
    a.set_ylabel("Diversity retained (%)")
    a.legend(fontsize=5.4)
    fs.panel_letter(a, "a")

    offs = {"Set1": -0.2, "Set2": 0.2}
    nmax = max(len(pd.read_csv(TAB / f"richness_{p.lower()}.csv")) for p in PANELS)
    for p in PANELS:
        r = pd.read_csv(TAB / f"richness_{p.lower()}.csv")
        x = np.arange(len(r)) + offs[p]
        b.bar(x, r["n_private_alleles"], 0.38, color=fs.PANEL_COL[p])
        c.bar(x, r["mean_rarefied_richness_Ag"], 0.38, color=fs.PANEL_COL[p])
    # single shared legend lives in panel (a) — b/c legends removed 2026-08-08 (the panel-b
    # copy overlapped the C5 bar; all three panels share the same two panel colours)
    for ax, ylab in ((b, "Private alleles (log scale)"),
                     (c, "Rarefied allelic richness (A$_g$)")):
        ax.set_xticks(np.arange(nmax))
        ax.set_xticklabels([fs.cluster_label(i) for i in range(nmax)], fontsize=5.4)
        ax.set_xlabel("Admixture cluster (per panel)")
        ax.set_ylabel(ylab)
    b.set_yscale("log")
    c.set_ylim(1.0, 2.0)
    fs.panel_letter(b, "b")
    fs.panel_letter(c, "c")

    fig.tight_layout()
    fs.save(fig, str(FIG_MAIN / "Fig08_application.png"))
    plt.close(fig)
    print("[35] Fig08_application written (3 panels)")


# ============================================================================
# Fig 9 — mechanism made visible (6 panels; needs stage 36 tables)
# ============================================================================
def fig09_mechanism_visible():
    need = ["height_by_cluster_set1.csv", "height_by_cluster_set2.csv",
            "height_by_cluster_summary.csv", "distance_pairs_set1.csv",
            "distance_pairs_set2.csv", "external_confusion_set1.csv",
            "classifier_baselines_set1.csv", "classifier_baselines_set2.csv",
            "external_label_classifier_set1.csv"]
    if not have(*need):
        print("[35] Fig09 skipped (stage 36 pending)")
        return
    fig, axes = plt.subplots(2, 3, figsize=(fs.DOUBLE, fs.DOUBLE * 0.66))
    (a, b, c), (d, e, f_) = axes
    kw = pd.read_csv(TAB / "height_by_cluster_summary.csv").set_index("panel")

    # (a,b) height index by admixture cluster — the confound, seen directly
    for ax, p, letter in ((a, "Set1", "a"), (b, "Set2", "b")):
        h = pd.read_csv(TAB / f"height_by_cluster_{p.lower()}.csv")
        clusters = sorted(h["cluster"].unique())
        data = [h.loc[h.cluster == cl, "height_index_z"] for cl in clusters]
        bp = ax.boxplot(data, positions=range(len(clusters)), widths=0.62,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color="0.15", lw=1.1),
                        whiskerprops=dict(lw=0.8), capprops=dict(lw=0.8))
        for patch, cl in zip(bp["boxes"], clusters):
            patch.set_facecolor(fs.cluster_color(int(cl)))
            patch.set_alpha(0.75)
            patch.set_edgecolor("0.25")
        rng = np.random.default_rng(7)
        for xi, vals in enumerate(data):
            ax.scatter(np.full(len(vals), xi) + rng.uniform(-0.16, 0.16, len(vals)),
                       vals, s=3.5, color="0.25", zorder=3, linewidths=0)
        ax.set_xticks(range(len(clusters)))
        ax.set_xticklabels([fs.cluster_label(int(cl)) for cl in clusters], fontsize=5.4)
        ax.set_xlabel("Admixture cluster")
        ax.set_ylabel("Height index (mean z-score)")
        kwp = kw.loc[p, "kruskal_p"]
        # short title — the long "Kruskal-Wallis" form collided with the next panel letter
        ax.set_title(f"{fs.PANEL_LABEL[p]} — KW P = {kwp:.0e}",
                     loc="left", fontsize=fs.FS_TITLE - 1)
        fs.panel_letter(ax, letter)

    # (c,d) the raw distance-pair cloud behind the primary Mantel r
    for ax, p, letter in ((c, "Set1", "c"), (d, "Set2", "d")):
        pairs = pd.read_csv(TAB / f"distance_pairs_{p.lower()}.csv")
        hb = ax.hexbin(pairs["genomic_dist"], pairs["phenomic_dist"], gridsize=38,
                       cmap="Blues", mincnt=1, linewidths=0.1)
        m = pd.read_csv(TAB / f"concordance_mantel_{p.lower()}.csv")
        r = float(m.loc[m.comparison == "genomic~phenomic", "r"].iloc[0])
        bt = pd.read_csv(TAB / f"mantel_bootstrap_{p.lower()}.csv").iloc[0]
        ax.annotate(f"Mantel r = {r:.3f}\n95% CI [{bt.ci95_lo:.3f}, {bt.ci95_hi:.3f}]",
                    xy=(0.03, 0.97), xycoords="axes fraction", va="top",
                    fontsize=fs.FS_ANNOT + 0.4)
        ax.set_xlabel("Genomic distance (1 − IBS)")
        ax.set_ylabel("Phenomic distance (Euclidean)")
        ax.set_title(fs.PANEL_LABEL[p], loc="left", fontsize=fs.FS_TITLE - 1)
        fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.04,
                     label="accession pairs").ax.tick_params(labelsize=fs.FS_TICK - 1)
        fs.panel_letter(ax, letter)

    # (e) external-label confusion matrix (Set 1)
    cm = pd.read_csv(TAB / "external_confusion_set1.csv")
    names = list(pd.unique(cm["true_group"]))
    M = cm.pivot(index="true_group", columns="pred_group",
                 values="count").reindex(index=names, columns=names).to_numpy()
    im = e.imshow(M, cmap="Blues")
    for i in range(len(names)):
        for j in range(len(names)):
            if M[i, j] > 0:
                e.text(j, i, int(M[i, j]), ha="center", va="center", fontsize=5.6,
                       color="white" if M[i, j] > 0.6 * M.max() else "0.1")
    e.set_xticks(range(len(names))); e.set_xticklabels(names, fontsize=5.4, rotation=45)
    e.set_yticks(range(len(names))); e.set_yticklabels(names, fontsize=5.4)
    e.set_xlabel("Predicted (from imaging alone)")
    e.set_ylabel("True 3K-RGP subpopulation")
    ext = pd.read_csv(TAB / "external_label_classifier_set1.csv").iloc[0]
    e.set_title(f"External labels — acc {ext.accuracy:.1%} "
                f"(null {ext.null_mean:.1%}, P = {ext.empirical_p:.3f})",
                loc="left", fontsize=fs.FS_TITLE - 1)
    fs.panel_letter(e, "e")

    # (f) robustness ladder: accuracy across target definitions
    ladder = [("Admixture (elbow K)", "admixture_elbowK", "accuracy"),
              ("k-means, K − 1", "kmeans_lo", "accuracy"),
              ("k-means, K + 1", "kmeans_hi", "accuracy"),
              ("Soft-Q weighted", "admixture_elbowK", "softQ_weighted_accuracy")]
    ys = np.arange(len(ladder) + 1)[::-1]
    for p in PANELS:
        cb = pd.read_csv(TAB / f"classifier_baselines_{p.lower()}.csv")
        full = cb[(cb.features == "full_204")]
        km = sorted(int(t.replace("kmeans_K", ""))
                    for t in full.target[full.target.str.startswith("kmeans")])
        vals = []
        for _, key, col in ladder:
            if key == "kmeans_lo":
                vals.append(float(full.loc[full.target == f"kmeans_K{km[0]}",
                                           "accuracy"].iloc[0]))
            elif key == "kmeans_hi":
                vals.append(float(full.loc[full.target == f"kmeans_K{km[1]}",
                                           "accuracy"].iloc[0]))
            else:
                vals.append(float(full.loc[full.target == "admixture_elbowK",
                                           col].iloc[0]))
        f_.plot(vals, ys[:len(ladder)], marker=fs.PANEL_MARKER[p], ms=5, lw=1.0,
                color=fs.PANEL_COL[p], label=fs.PANEL_LABEL[p])
        null = float(full.loc[full.target == "admixture_elbowK", "null_mean"].iloc[0])
        f_.axvline(null, lw=0.8, ls=":", color=fs.PANEL_COL[p], alpha=0.8)
    ext = pd.read_csv(TAB / "external_label_classifier_set1.csv").iloc[0]
    f_.scatter([ext.accuracy], [ys[len(ladder)]], marker="D", s=30, color="#0f5f59",
               zorder=5)
    f_.scatter([ext.null_mean], [ys[len(ladder)]], marker="D", s=30,
               facecolors="none", edgecolors="#0f5f59", zorder=5)
    f_.set_yticks(ys)
    f_.set_yticklabels([lab for lab, _, _ in ladder] + ["External 3K-RGP\n(Set 1; open = null)"],
                       fontsize=5.4)
    f_.set_xlabel("Classification accuracy (dotted = null)")
    f_.legend(fontsize=5.2, loc="lower left")   # lower right collided with the null diamond
    fs.panel_letter(f_, "f")

    fig.tight_layout()
    fs.save(fig, str(FIG_MAIN / "Fig09_mechanism_visible.png"))
    plt.close(fig)
    print("[35] Fig09_mechanism_visible written (6 panels)")


# ============================================================================
# Fig 10 — two platforms, one biology (3 panels; needs stage 36 scoreboard)
# ============================================================================
def fig10_platforms():
    need = ["replication_scoreboard.csv", "maf_spectrum_set1.csv", "maf_spectrum_set2.csv",
            "snp_density_by_chrom.csv"]
    if not have(*need):
        print("[35] Fig10 skipped (stage 36 pending)")
        return
    fig, axes = plt.subplots(1, 3, figsize=(fs.DOUBLE, fs.DOUBLE * 0.40),
                             gridspec_kw={"width_ratios": [1.3, 1, 1], "wspace": 0.45})
    a, b, c = axes

    # (a) replication scoreboard dumbbell
    sc = pd.read_csv(TAB / "replication_scoreboard.csv")
    ys = np.arange(len(sc))[::-1]
    for yi, (_, row) in zip(ys, sc.iterrows()):
        a.plot([row.set1, row.set2], [yi, yi], color="0.75", lw=1.0, zorder=1)
    a.scatter(sc.set1, ys, s=22, marker=fs.PANEL_MARKER["Set1"],
              color=fs.PANEL_COL["Set1"], zorder=3, label=fs.PANEL_LABEL["Set1"])
    a.scatter(sc.set2, ys, s=22, marker=fs.PANEL_MARKER["Set2"],
              color=fs.PANEL_COL["Set2"], zorder=3, label=fs.PANEL_LABEL["Set2"])
    a.set_yticks(ys)
    a.set_yticklabels(sc.metric, fontsize=5.2)
    a.set_xlabel("Estimate (dimensionless)")
    # centre-right rows (Gst/theta/Mantel) top out at ~0.45 — verified empty; lower-right sat
    # between the ARI and core-retention markers (QA sweep)
    a.legend(fontsize=5.2, loc="center right")
    fs.panel_letter(a, "a")

    # (b) MAF spectrum: WGS vs ascertained array
    w = 0.38
    m1 = pd.read_csv(TAB / "maf_spectrum_set1.csv")
    bins = list(m1["maf_bin"])
    x = np.arange(len(bins))
    for i, p in enumerate(PANELS):
        mm = pd.read_csv(TAB / f"maf_spectrum_{p.lower()}.csv").set_index("maf_bin")
        share = (mm["n_snps"] / mm["n_snps"].sum()).reindex(bins)
        b.bar(x + (i - 0.5) * w, share, width=w, color=fs.PANEL_COL[p])
    b.set_xticks(x)
    b.set_xticklabels(bins, fontsize=5.2, rotation=45)
    b.set_xlabel("Minor-allele-frequency bin (QC'd SNPs)")
    b.set_ylabel("Share of SNPs")
    fs.panel_letter(b, "b")

    # (c) SNP density per chromosome (log): why array replication is nontrivial
    dd = pd.read_csv(TAB / "snp_density_by_chrom.csv")
    for i, p in enumerate(PANELS):
        sub = dd[dd.panel == p].sort_values("chrom")
        c.bar(sub["chrom"] + (i - 0.5) * 0.38, sub["snps_per_mb"], width=0.38,
              color=fs.PANEL_COL[p])
    c.set_yscale("log")
    c.set_xticks(range(1, 13))
    c.set_xlabel("Chromosome")
    c.set_ylabel("QC'd SNPs per Mb (log scale)")
    fs.panel_letter(c, "c")

    fig.tight_layout()
    fs.save(fig, str(FIG_MAIN / "Fig10_platforms.png"))
    plt.close(fig)
    print("[35] Fig10_platforms written (3 panels)")


if __name__ == "__main__":
    fig02_structure()
    fig04_diff_diversity_ld()
    fig05_concordance_mechanism()
    fig08_application()
    fig09_mechanism_visible()
    fig10_platforms()
    print("[35] done")
