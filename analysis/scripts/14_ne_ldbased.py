#!/usr/bin/env python
"""14_ne_ldbased.py — [stage 14] linkage-disequilibrium-based effective population size (Ne),
with an explicit selfing adjustment and an explicit assumption ledger.

Motivation: Ne is the policy-relevant number behind
the narrow-base/core-collection message. Stage 06 deliberately refused to hand-roll Ne without a
verifiable implementation; this stage closes that gap by implementing the textbook estimator with
every assumption stated and a sensitivity grid, and by presenting the result as an
order-of-magnitude trajectory, not a precise point estimate.

Method (Sved 1971; sample-size correction as in Weir & Hill 1980; binned-distance trajectory as
popularised by SNeP, Corbin et al. 2012):
    E[r^2] ~= 1/(1 + 4*Ne*c) + 1/n
 => Ne(c)  = (1/(4c)) * (1/(r2_bin - 1/n) - 1),  at time t ~= 1/(2c) generations ago,
using the stage-07 binned mean r^2 (25-kb bins to 1 Mb). Physical -> genetic distance uses the
genome-wide average rice map density (~4.0 cM/Mb; Harushima et al. 1998 map length ~1,521 cM over
~373 Mb), with a 3-5 cM/Mb sensitivity band.

Selfing adjustment: partial self-fertilisation reduces EFFECTIVE recombination,
c* = c * (1 - F) (Nordborg 2000 scaling; F from stage 16). Under near-complete selfing
(F ~ 0.93/0.79) the panmictic formula therefore UNDERSTATES Ne by ~1/(1-F); we report both the
panmictic and the selfing-adjusted trajectories and treat the mating system as the dominant
uncertainty. Caveats reported with the result: equilibrium assumptions, admixed panel (structure
inflates LD and biases Ne downward), MAF>0.05 ascertainment of the r^2 input, Set1 r^2 computed
on a 10%-thinned marker set.

Outputs (analysis/results/tables/):
  ne_trajectory_{set}.csv — bin_mid_bp, c_morgan, t_gen, r2, r2_adj, ne_panmictic, ne_selfing_adj
  ne_summary.csv          — per panel: recent-Ne summary (harmonic mean over the 3 most-recent
                             usable bins), both adjustments, cM/Mb sensitivity, assumption string
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import TAB

CM_PER_MB = 4.0                 # genome-wide average; sensitivity band below
CM_SENS = (3.0, 5.0)
ALPHA = 1.0                     # Sved drift-only; mutation-corrected alpha=2.2 noted in Methods
N_IND = {"Set1": 150, "Set2": 147}
MIN_BIN_BP = 25_000             # drop the shortest bin (gene conversion / strong short-range LD)


def trajectory(panel, cm_per_mb):
    ld = pd.read_csv(TAB / f"ld_decay_{panel.lower()}.csv")
    ld = ld[(ld["panel"] == panel) & (ld["bin_mid_bp"] >= MIN_BIN_BP)].copy()
    n = N_IND[panel]
    F = pd.read_csv(TAB / "hwe_summary.csv").set_index("panel").loc[panel, "mean_F"]
    c = ld["bin_mid_bp"].to_numpy() * cm_per_mb * 1e-8          # Morgans
    r2 = ld["mean_r2"].to_numpy()
    r2_adj = r2 - 1.0 / n
    ok = r2_adj > 0.01                                          # avoid exploding 1/r2_adj
    ne_pan = np.where(ok, (1.0 / (4 * c)) * (1.0 / r2_adj - ALPHA), np.nan)
    c_star = c * (1.0 - F)
    ne_self = np.where(ok, (1.0 / (4 * c_star)) * (1.0 / r2_adj - ALPHA), np.nan)
    out = pd.DataFrame({
        "bin_mid_bp": ld["bin_mid_bp"], "n_pairs": ld["n_pairs"],
        "c_morgan": c, "t_gen_panmictic": 1.0 / (2 * c), "t_gen_selfing_adj": 1.0 / (2 * c_star),
        "r2": r2, "r2_adj": r2_adj,
        "ne_panmictic": ne_pan, "ne_selfing_adj": ne_self})
    return out, float(F)


def recent_ne(traj, col):
    # harmonic mean over the 3 largest-distance (most recent time) usable bins
    v = traj.dropna(subset=[col]).sort_values("bin_mid_bp", ascending=False)[col].head(3)
    return float(len(v) / (1.0 / v).sum()) if len(v) else np.nan


def main():
    rows = []
    for panel in ["Set1", "Set2"]:
        traj, F = trajectory(panel, CM_PER_MB)
        traj.to_csv(TAB / f"ne_trajectory_{panel.lower()}.csv", index=False)
        row = {"panel": panel, "F_used": F,
               "recent_ne_panmictic": recent_ne(traj, "ne_panmictic"),
               "recent_ne_selfing_adj": recent_ne(traj, "ne_selfing_adj")}
        for cm in CM_SENS:
            t2, _ = trajectory(panel, cm)
            row[f"recent_ne_panmictic_cM{cm:g}"] = recent_ne(t2, "ne_panmictic")
        row["assumptions"] = (f"Sved1971; alpha={ALPHA}; r2_adj=r2-1/n; {CM_PER_MB} cM/Mb; "
                              f"selfing c*=c(1-F); bins>={MIN_BIN_BP}bp; "
                              "structure+MAF ascertainment biases stated in Methods")
        rows.append(row)
        print(f"[14] {panel}: recent Ne (panmictic) ~ {row['recent_ne_panmictic']:.0f}; "
              f"selfing-adjusted ~ {row['recent_ne_selfing_adj']:.0f} (F = {F:.3f})")
    pd.DataFrame(rows).to_csv(TAB / "ne_summary.csv", index=False)
    print("[14] done -> ne_summary.csv")


if __name__ == "__main__":
    main()
