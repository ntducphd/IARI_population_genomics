#!/usr/bin/env python
"""figstyle.py — shared visual system for the compendium figures (house convention:
journal column widths in mm, luminance-spaced categorical palettes so nothing is
colour-alone, PNG for review + vector PDF for the journal from one call).
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

MM = 1.0 / 25.4
SINGLE = 85 * MM
ONE_HALF = 130 * MM
DOUBLE = 180 * MM

PANEL_COL = {"Set1": "#2b4a7d", "Set2": "#c25f2a"}     # WGS vs 50K array, spaced on luminance
PANEL_LABEL = {"Set1": "Set 1 (WGS)", "Set2": "Set 2 (50K array)"}
PANEL_MARKER = {"Set1": "o", "Set2": "^"}

# Admixture-cluster palette, shared by every cluster-coloured panel (Fig 2 PCA, Fig 3 Q bars,
# Fig 5 heatmap labels, Fig 9b): Okabe-Ito colour-blind-safe set (8) + a 9th neutral, indexed
# C0..C8 to MATCH the C-labels used in the text and every table (never "cluster 1..9").
CLUSTER_COL = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
               "#D55E00", "#CC79A7", "#999999", "#000000"]


def cluster_color(i):
    return CLUSTER_COL[i % len(CLUSTER_COL)]


def cluster_label(i):
    return f"C{i}"

TREAT = {"Control": "#0f5f59", "NStress": "#e8912e"}
FEATGROUP_COL = {"Colour": "#9b8ec4", "NIR": "#c25f2a", "Size_morphology": "#2b4a7d"}

FS_TITLE, FS_LABEL, FS_TICK = 8.0, 7.5, 6.5
FS_ANNOT, FS_LEGEND, FS_PANEL = 6.5, 6.5, 9.0


def set_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
        "font.size": 7.5, "axes.titlesize": FS_TITLE, "axes.titleweight": "bold",
        "axes.labelsize": FS_LABEL, "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND, "axes.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.2, "ytick.major.size": 2.2,
        "lines.linewidth": 1.0, "patch.linewidth": 0.5,
        "legend.frameon": False, "legend.handlelength": 1.4, "legend.columnspacing": 1.0,
        "figure.dpi": 300, "savefig.dpi": 600, "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "mathtext.default": "regular",
    })


def panel_letter(ax, letter, dx=-24, dy=5, size=None):
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction",
                xytext=(dx, dy), textcoords="offset points",
                fontsize=size or FS_PANEL, fontweight="bold", va="bottom", ha="left")


def star(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."


def save(fig, path):
    """PNG for review and PDF for the journal, from one call. bbox_inches='tight' ensures anything
    drawn slightly outside the nominal figure bounds (e.g. a suptitle at y>1.0) is not clipped --
    found clipping the top of several titles before this fix (Fig1/3/4/7, SuppFig1/2)."""
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(path.rsplit(".", 1)[0] + ".pdf", bbox_inches="tight")
    # SVG with real text elements (svg.fonttype = "none"): the editable-figure source for
    # the journal figure booklets — Word can convert these to native editable shapes/text
    fig.savefig(path.rsplit(".", 1)[0] + ".svg", bbox_inches="tight")