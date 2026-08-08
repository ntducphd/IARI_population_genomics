#!/usr/bin/env python
"""25_ld_check.py — [stage 25] LD-decay reporting upgrades .

Two fixes to the stage-07 LD analysis:
  [A] Standard-threshold decay distances alongside the bespoke half-of-shortest-bin definition:
      the distance at which the binned mean r^2 curve first crosses r^2 = 0.2 and 0.1 (linear
      interpolation between bin midpoints) — the conventions most rice LD papers report — so the
      manuscript's 562 kb / 212 kb numbers become directly comparable with the cited literature.
  [B] An EMPIRICAL thinning-bias check: stage 07 thinned Set 1 to ~10% for tractability and the
      manuscript claimed (uncited) that thinning "likely modestly over-estimates" half-decay.
      Here chromosome 1, window 1-10 Mb, is computed UNTHINNED and again 10%-thinned (same seed
      protocol), and both curves' half-decay and r2=0.2 crossings are compared — replacing the
      hand-waved direction-of-bias claim with a measured one.

Outputs (analysis/results/tables/):
  ld_threshold_summary.csv — per panel: bespoke half-decay (stage 07), dist at r2=0.2, r2=0.1
  ld_thinning_check.csv    — chr1 1-10 Mb: unthinned vs thinned half-decay + r2 crossings
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PLINK, INTERIM, TAB

BIN_KB = 25


def run_plink(*args):
    cmd = [str(PLINK), *[str(a) for a in args], "--chr-set", "12", "no-xy",
           "--allow-extra-chr", "--silent"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("PLINK ERROR:\n", r.stderr[-1500:])
        raise SystemExit(1)


def crossing(bins_bp, r2, threshold):
    """First distance where the binned curve crosses `threshold` (linear interpolation)."""
    for i in range(1, len(r2)):
        if r2[i - 1] >= threshold > r2[i]:
            x0, x1, y0, y1 = bins_bp[i - 1], bins_bp[i], r2[i - 1], r2[i]
            return float(x0 + (y0 - threshold) * (x1 - x0) / (y0 - y1))
    return np.nan


def curve_stats(df):
    df = df.sort_values("bin_mid_bp")
    b, r = df["bin_mid_bp"].to_numpy(), df["mean_r2"].to_numpy()
    half = crossing(b, r, r[0] / 2)
    return {"r2_shortest_bin": float(r[0]), "half_decay_bp": half,
            "dist_r2_0.2_bp": crossing(b, r, 0.2), "dist_r2_0.1_bp": crossing(b, r, 0.1)}


def binned_from_ldfile(ld_path):
    ld = pd.read_csv(ld_path, sep=r"\s+")
    d = (ld["BP_B"] - ld["BP_A"]).abs()
    bins = (d // (BIN_KB * 1000)).astype(int)
    g = ld.groupby(bins)["R2"].agg(["mean", "size"]).reset_index()
    g["bin_mid_bp"] = g["index"] * BIN_KB * 1000 + BIN_KB * 500
    return pd.DataFrame({"bin_mid_bp": g["bin_mid_bp"], "mean_r2": g["mean"],
                         "n_pairs": g["size"]})


def main():
    # [A] thresholds on the existing stage-07 binned curves
    rows = []
    for panel in ["Set1", "Set2"]:
        cur = pd.read_csv(TAB / f"ld_decay_{panel.lower()}.csv")
        st = curve_stats(cur[cur["panel"] == panel])
        st["panel"] = panel
        rows.append(st)
        print(f"[25A] {panel}: half-decay {st['half_decay_bp']/1e3:.0f} kb; "
              f"r2=0.2 at {st['dist_r2_0.2_bp']/1e3:.0f} kb; "
              f"r2=0.1 at {st['dist_r2_0.1_bp']/1e3:.0f} kb")
    pd.DataFrame(rows).to_csv(TAB / "ld_threshold_summary.csv", index=False)

    # [B] chr1 1-10 Mb, unthinned vs 10%-thinned
    checks = []
    for tag, extra in [("unthinned", []), ("thinned10pct", ["--thin", "0.1", "--seed", "42"])]:
        out = INTERIM / f"ldchk_{tag}"
        run_plink("--bfile", INTERIM / "Set1_qc", "--chr", "1", "--from-bp", "1",
                  "--to-bp", "10000000", *extra, "--r2", "--ld-window-kb", "1000",
                  "--ld-window", "99999", "--ld-window-r2", "0", "--out", out)
        cur = binned_from_ldfile(f"{out}.ld")
        st = curve_stats(cur)
        st["variant"] = tag
        checks.append(st)
        print(f"[25B] chr1 1-10Mb {tag}: half-decay {st['half_decay_bp']/1e3:.0f} kb, "
              f"r2=0.2 at {st['dist_r2_0.2_bp']/1e3:.0f} kb")
    pd.DataFrame(checks).to_csv(TAB / "ld_thinning_check.csv", index=False)
    print("[25] done")


if __name__ == "__main__":
    main()
