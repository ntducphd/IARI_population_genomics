#!/usr/bin/env python
"""29_fig01_design.py — Figure 1 built to the house schematic standard:
numbered stage badges, rounded boxes with bold head + one detail line, semantic colour coding
with a legend, converging non-crossing arrows, uniform type hierarchy, real numbers inline.

Balance pass 2026-08-08 (v3): one 3-column grid shared by rows 1-2 (x = 9/45.5/82, w = 33), an
equal 4-column grid for row 3 (w = 24.25), every vertical arrow lands on its TARGET box centre,
the Set1->anchoring link runs as an elbow in the inter-row gap with a hop where it crosses the
Set2 arrow, Pillar B in/out arrows are symmetric about the canvas centre (x = 62), and the
vertical rhythm is a constant 4.6-unit gap between rows.

Output: figures/main/Fig01_design.png + .pdf (the bundle renumbers to Figure1).
"""
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import FIG_MAIN
import figstyle as fs

fs.set_style()
mpl.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})

FS, FS_ST, FS_SUP, FS_LEG = 8.6, 10.5, 13.0, 8.6
S1C = fs.PANEL_COL["Set1"]          # WGS blue
S2C = fs.PANEL_COL["Set2"]          # 50K orange
PHC = "#0f5f59"                      # phenomics teal
NEUT, EDGE, INK = "#F2F5F9", "#33404F", "#161616"
BS = "round,pad=0.02,rounding_size=0.9"

# shared grids
C3X, C3W = (9, 45.5, 82), 33                      # rows 1-2: three columns
C4X, C4W = (9, 35.5, 62, 88.5), 24.25             # row 3: four columns
C3C = [x + C3W / 2 for x in C3X]                  # 25.5, 62, 98.5
C4C = [x + C4W / 2 for x in C4X]                  # 21.125, 47.625, 74.125, 100.625


def tint(c, f=0.12):
    r, g, b = to_rgb(c)
    return (r * f + (1 - f), g * f + (1 - f), b * f + (1 - f))


fig, ax = plt.subplots(figsize=(12.6, 7.9))
ax.set_xlim(0, 124)
ax.set_ylim(7.5, 74.5)
ax.axis("off")


def box(x, y, w, h, title, detail="", fc=NEUT, ec=EDGE, tc=INK, lw=1.4):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=BS, fc=fc, ec=ec, lw=lw, zorder=3))
    if detail:
        ax.text(x + w / 2, y + h * 0.64, title, ha="center", va="center",
                fontsize=FS, fontweight="bold", color=tc, zorder=4)
        ax.text(x + w / 2, y + h * 0.28, detail, ha="center", va="center",
                fontsize=FS, color="#2A2A2A", zorder=4)
    else:
        ax.text(x + w / 2, y + h / 2, title, ha="center", va="center",
                fontsize=FS, fontweight="bold", color=tc, zorder=4)


def arrow(x1, y1, x2, y2, c=EDGE, lw=1.6, sb=1.6):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14,
                                 lw=lw, color=c, zorder=2, shrinkA=1.2, shrinkB=sb))


def stage(n, y, t, cx=3.2, tx=6.0):
    ax.add_patch(Circle((cx, y), 1.45, fc=EDGE, ec="none", zorder=5))
    ax.text(cx, y, str(n), ha="center", va="center", fontsize=8.4, color="white",
            fontweight="bold", zorder=6)
    ax.text(tx, y, t, ha="left", va="center", fontsize=FS_ST, fontweight="bold",
            color=EDGE, zorder=5)


# ---- legend (no in-image title — the figure number/title belongs to the legend text,
# QA pass 2026-08-08) ---------------------------------------------------------
for i, (lab, c) in enumerate([("Set 1 (WGS)", S1C), ("Set 2 (50K array)", S2C),
                              ("Phenomics", PHC)]):
    x0 = 74 + i * 16
    ax.add_patch(FancyBboxPatch((x0, 71.4), 2.6, 1.7, boxstyle=BS, fc=tint(c), ec=c,
                                lw=1.2, zorder=5))
    ax.text(x0 + 3.4, 72.25, lab, ha="left", va="center", fontsize=FS_LEG, zorder=5)

# ---- stage 1: panels --------------------------------------------------------
stage(1, 72.3, "Experimental data")
box(9, 65.2, 106, 4.9, "Two disjoint rice panels — same breeding programme",
    "no shared accessions · single greenhouse season · canopy-imaged and phenotyped")
box(C3X[0], 55.6, C3W, 6.8, "Set 1 — WGS",
    "n = 150 · 1.01 M raw → 502,675 QC → 24,370 pruned", ec=S1C, fc=tint(S1C, 0.07))
box(C3X[1], 55.6, C3W, 6.8, "Set 2 — 50K array",
    "n = 147 · 50,051 raw → 31,565 QC → 2,942 pruned", ec=S2C, fc=tint(S2C, 0.07))
box(C3X[2], 55.6, C3W, 6.8, "Canopy imaging + traits",
    "204 features × 3 timepoints · agronomic / NUE traits", ec=PHC, fc=tint(PHC, 0.07))
for xc in C3C:
    arrow(xc, 65.2, xc, 62.4)

# ---- stage 2: structure & ancestry (same 3-column grid as row 1) ------------
stage(2, 53.0, "Structure & ancestry")
box(C3X[0], 44.2, C3W, 6.8, "Four concordant methods",
    "PCA · sNMF · UMAP · DAPC — ARI 0.60–0.76")
box(C3X[1], 44.2, C3W, 6.8, "Trees & kinship", "NJ (bootstrap 100×) · GRM")
box(C3X[2], 44.2, C3W, 6.8, "Global anchoring", "3K-RGP labels · f3 admixture tests",
    ec=S1C, fc=tint(S1C, 0.05))
arrow(C3C[0], 55.6, C3C[0], 51.0, c=S1C)          # Set 1 -> structure (centre to centre)
arrow(C3C[1], 55.6, C3C[1], 51.0, c=S2C)          # Set 2 -> trees/kinship
# Set 1 -> global anchoring: elbow in the inter-row gap, hop where it crosses the Set2 arrow
YB = 53.3
ax.plot([38, 38, 60.9], [55.6, YB, YB], color=S1C, lw=1.6, zorder=2,
        solid_capstyle="round", solid_joinstyle="round")
ax.plot([63.1, C3C[2]], [YB, YB], color=S1C, lw=1.6, zorder=2,
        solid_capstyle="round", solid_joinstyle="round")
arrow(C3C[2], YB, C3C[2], 51.0, c=S1C)

# ---- stage 3: diversity, history, selection (equal 4-column grid) -----------
stage(3, 41.3, "Genome layers")
box(C4X[0], 32.8, C4W, 6.8, "Diversity & mating", "He/Ho · ŝ = 0.88–0.96 · richness")
box(C4X[1], 32.8, C4W, 6.8, "Genome history", "ROH/F_ROH · LD-Ne · Stairway")
box(C4X[2], 32.8, C4W, 6.8, "Differentiation", "WC θ (CI) · AMOVA ×2 · LD decay")
box(C4X[3], 32.8, C4W, 6.8, "Selection scans", "pcadapt · iHS / XP-EHH")
for xc in C4C:
    arrow(xc, 44.2, xc, 39.6)

# ---- stage 4: pillar B (punch line) ----------------------------------------
stage(4, 29.8, "Concordance — Pillar B")
box(24, 21.0, 76, 7.2, "Does canopy imaging carry genetic ancestry?",
    "Mantel / partial Mantel · Procrustes · RF vs permutation null · MMRR — "
    "replicated in both disjoint panels",
    ec=PHC, fc=tint(PHC, 0.14), lw=2.0)
arrow(C4C[1], 32.8, 54, 28.2)                     # symmetric in-arrows about x = 62
arrow(C4C[2], 32.8, 70, 28.2)
# phenomics feeds Pillar B — routed along the right margin, clear of every box
ax.plot([115, 121, 121], [59, 59, 24.6], color=PHC, lw=1.6, zorder=2,
        solid_capstyle="round", solid_joinstyle="round")
arrow(121, 24.6, 100.8, 24.6, c=PHC)

# ---- stage 5: outputs -------------------------------------------------------
stage(5, 18.0, "Outputs")
box(9, 9.2, 50, 7.2, "Population-structure reference",
    "stratification corrections + disjoint-panel logic for the programme's studies")
box(64, 9.2, 51, 7.2, "Breeding resource",
    "core collection (10%) · Ne / F_ROH management · underused accessions")
arrow(48, 21.0, 34, 16.4)                         # symmetric out-arrows, land on box centres
arrow(76, 21.0, 89.5, 16.4)

fig.tight_layout()
fs.save(fig, str(FIG_MAIN / "Fig01_design.png"))
plt.close(fig)
print("[29] Fig01_design rebuilt (balanced grid, house schematic standard)")
