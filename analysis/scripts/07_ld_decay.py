#!/usr/bin/env python
"""07_ld_decay.py — LD decay (r^2 vs physical distance) per panel, half-decay distance.

Set1 (WGS, 502,675 QC SNPs) is thinned to a random ~10% subsample first (`plink --thin 0.1`) before
computing pairwise r^2 within a 1 Mb window: at full density, a 1 Mb window around every one of
502k SNPs would emit on the order of 10^8 pairs (average SNP spacing ~800 bp * ~1250 neighbours
each side), which is not a tractable file size/runtime on this machine and is not necessary --
thinning is standard practice for LD-decay curves (the decay pattern, not every single pair,
is what the curve needs) and keeps physical positions intact. Set2 (50K array, 2,942 QC SNPs) is
already sparse enough to use directly, no thinning.

Half-decay distance: the smallest bin midpoint at which mean r^2 first drops to <= half of the
mean r^2 in the shortest-distance bin (the same simple, standard definition used across the rice
LD-decay literature, e.g. the ~196 kb figure cited for DRC rice DArTseq panels in
the literature-derived competitor list).

Outputs (analysis/results/tables/):
  ld_decay_{set}.csv          — bin_mid_bp, mean_r2, n_pairs
  ld_decay_summary.csv        — one row per panel: half_decay_bp, r2_at_shortest_bin, n_pairs_total
"""
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PLINK, INTERIM, TAB

WINDOW_KB = 1000
BIN_KB = 25
PANELS = {"Set1": True, "Set2": False}   # value = whether to thin first


def run_plink(*args):
    cmd = [str(PLINK), *[str(a) for a in args], "--chr-set", "12", "no-xy",
           "--allow-extra-chr", "--silent"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("PLINK ERROR:\n", r.stderr[-2000:])
        raise SystemExit(1)
    return r


rows_summary = []
for panel, do_thin in PANELS.items():
    print(f"\n=== [{panel}] LD decay ===")
    stem = INTERIM / f"{panel}_qc"
    ldir = INTERIM / "lea"; ldir.mkdir(parents=True, exist_ok=True)

    if do_thin:
        thin_stem = ldir / f"{panel}_ld_thin"
        run_plink("--bfile", stem, "--thin", "0.1", "--make-bed", "--out", thin_stem, "--seed", "42")
        ld_input = thin_stem
    else:
        ld_input = stem

    ld_out = ldir / f"{panel}_ld"
    run_plink("--bfile", ld_input, "--r2", "--ld-window-kb", str(WINDOW_KB),
              "--ld-window", "999999", "--ld-window-r2", "0", "--out", ld_out)

    ld_file = ld_out.with_suffix(".ld")
    n_input_snps = sum(1 for _ in open(f"{ld_input}.bim"))
    print(f"  input: {n_input_snps} SNPs -> {ld_file.name}")

    ld = pd.read_csv(ld_file, sep=r"\s+")
    ld["dist_bp"] = (ld["BP_B"] - ld["BP_A"]).abs()
    ld = ld[ld["dist_bp"] > 0]
    print(f"  {len(ld)} SNP pairs within {WINDOW_KB} kb")

    bin_edges = np.arange(0, WINDOW_KB * 1000 + BIN_KB * 1000, BIN_KB * 1000)
    ld["bin"] = pd.cut(ld["dist_bp"], bins=bin_edges, right=False)
    binned = ld.groupby("bin", observed=True).agg(mean_r2=("R2", "mean"), n_pairs=("R2", "size")).reset_index()
    binned["bin_mid_bp"] = binned["bin"].apply(lambda b: int(b.mid))
    binned = binned[["bin_mid_bp", "mean_r2", "n_pairs"]].sort_values("bin_mid_bp")
    binned.insert(0, "panel", panel)
    binned.to_csv(TAB / f"ld_decay_{panel.lower()}.csv", index=False)

    r2_shortest = binned.iloc[0]["mean_r2"]
    half = r2_shortest / 2.0
    below_half = binned[binned["mean_r2"] <= half]
    half_decay_bp = int(below_half.iloc[0]["bin_mid_bp"]) if len(below_half) else None
    print(f"  r2 at shortest bin ({binned.iloc[0]['bin_mid_bp']} bp) = {r2_shortest:.4f}")
    if half_decay_bp is not None:
        print(f"  half-decay distance ~ {half_decay_bp:,} bp ({half_decay_bp/1000:.0f} kb)")
    else:
        print(f"  half-decay distance NOT reached within {WINDOW_KB} kb (r2 stays above half-max)")

    rows_summary.append(dict(panel=panel, thinned=do_thin, n_snps_used=n_input_snps,
                              n_pairs_total=int(binned["n_pairs"].sum()),
                              r2_at_shortest_bin=round(float(r2_shortest), 4),
                              half_decay_bp=half_decay_bp))

pd.DataFrame(rows_summary).to_csv(TAB / "ld_decay_summary.csv", index=False)
print("\n-> tables/ld_decay_{set}.csv + ld_decay_summary.csv")
