#!/usr/bin/env python
"""09b_core_collection.py — Stage 9, Pillar C breeding resource: a core/mini-core collection that
maximises genetic diversity coverage, plus a list of underused (diverse-tail) accessions.

Why Python, not corehunter (R)? corehunter installs cleanly from CRAN but its rJava backend
crashes R deterministically (`rJava::.jinit()` reproducibly reproduced an access violation on this
machine, retested -- see the environment-audit record §3). Per the agreed fallback, this
implements the classic **M (maximisation) strategy** directly: greedy max-min diversity / farthest-
point sampling on the genomic distance matrix (Schoen & Brown 1993; Franco et al. 2005; this is the
same family of algorithm corehunter itself implements under its "EN"/max-diversity objectives, just
without the Java dependency):
  1. seed the core with the most distant pair of accessions
  2. repeatedly add the accession whose MINIMUM distance to the current core is largest
     (this spreads the core across the full diversity space rather than clustering it)
  3. do this for a range of core sizes and report the diversity RETAINED at each size (mean
     pairwise distance of the core / mean pairwise distance of the full panel), the standard
     core-collection evaluation criterion
  4. recommend the smallest core retaining >= 80% of full-panel diversity (a conventional target)

Outputs (analysis/results/tables/):
  core_collection_{set}.csv          — accession, in_core (at the recommended size), core_rank
                                        (the order it was added -- rank 1/2 are the seed pair)
  core_collection_curve_{set}.csv     — core_size, pct_diversity_retained
  core_collection_summary_{set}.csv  — recommended_size, pct_of_panel, pct_diversity_retained
  underused_accessions_{set}.csv     — accessions NOT in the recommended core, sorted by how much
                                        they'd add if included next (a pre-breeding priority list)
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import TAB

CORE_SIZE_GRID_PCT = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]
PANELS = ["Set1", "Set2"]


def load_dist_matrix(panel):
    ibs = pd.read_csv(TAB / f"ibs_dist_{panel.lower()}.csv")
    ids = sorted(set(ibs["sample_a"]) | set(ibs["sample_b"]))
    p = ibs.pivot(index="sample_a", columns="sample_b", values="distance").reindex(index=ids, columns=ids)
    D = p.to_numpy().copy()
    np.fill_diagonal(D, 0.0)
    return ids, D


def maxmin_core_order(D, ids):
    """Greedy max-min diversity (farthest-point sampling): returns accessions in the order they
    would be added to the core, starting from the single most-distant pair."""
    n = len(ids)
    i, j = np.unravel_index(np.argmax(D), D.shape)
    order = [i, j]
    in_core = np.zeros(n, dtype=bool)
    in_core[[i, j]] = True
    min_dist_to_core = np.minimum(D[i], D[j])
    while len(order) < n:
        min_dist_to_core[in_core] = -np.inf
        nxt = int(np.argmax(min_dist_to_core))
        order.append(nxt)
        in_core[nxt] = True
        min_dist_to_core = np.minimum(min_dist_to_core, D[nxt])
    return order


def mean_upper(D, idx=None):
    sub = D if idx is None else D[np.ix_(idx, idx)]
    iu = np.triu_indices_from(sub, k=1)
    return float(sub[iu].mean()) if len(iu[0]) else float("nan")


for panel in PANELS:
    print(f"\n=== [{panel}] core collection (M-strategy / max-min diversity) ===")
    ids, D = load_dist_matrix(panel)
    n = len(ids)
    full_mean_dist = mean_upper(D)
    order = maxmin_core_order(D, ids)
    print(f"  n={n} accessions, full-panel mean pairwise distance = {full_mean_dist:.4f}")

    curve_rows = []
    for pct in CORE_SIZE_GRID_PCT:
        k = max(2, round(pct * n))
        idx = order[:k]
        retained = mean_upper(D, idx) / full_mean_dist if full_mean_dist > 0 else float("nan")
        curve_rows.append(dict(panel=panel, core_size=k, pct_of_panel=round(100 * k / n, 1),
                                pct_diversity_retained=round(100 * retained, 1)))
        print(f"    core size {k:3d} ({100*k/n:4.1f}% of panel): diversity retained = {100*retained:.1f}%")

    curve_df = pd.DataFrame(curve_rows)
    curve_df.to_csv(TAB / f"core_collection_curve_{panel.lower()}.csv", index=False)

    # Recommendation: the conventional ~10%-of-collection core size (Frankel & Brown 1984; van
    # Hintum et al. 2000 cite 5-20% as the typical range), NOT "smallest size clearing a diversity
    # threshold". Max-min/farthest-point sampling picks extremes first BY CONSTRUCTION, so mean
    # pairwise distance in the core exceeds the full-panel mean well above 100% at very small
    # sizes (confirmed above: 105-113% retained already at ~5% of the panel) -- a threshold-
    # crossing rule degenerates to "use almost no accessions" and is not a usable breeding
    # recommendation. The retention curve itself (written above) remains the transparent,
    # decision-relevant artefact; this just picks a practically sized point on it.
    target_pct = 10.0
    rec_idx = (curve_df["pct_of_panel"] - target_pct).abs().idxmin()
    rec = curve_df.loc[rec_idx]
    rec_size = int(rec["core_size"])
    print(f"  recommended core: {rec_size} accessions (~{target_pct:.0f}% of panel, the conventional "
          f"core-collection size per Frankel & Brown 1984 / van Hintum et al. 2000), "
          f"{rec['pct_diversity_retained']:.1f}% diversity retained (full curve in core_collection_curve_*.csv)")

    rank_of = {ids[o]: r + 1 for r, o in enumerate(order)}
    core_df = pd.DataFrame({"panel": panel, "accession": ids})
    core_df["in_core"] = [rank_of[a] <= rec_size for a in ids]
    core_df["core_rank"] = [rank_of[a] for a in ids]
    core_df = core_df.sort_values("core_rank")
    core_df.to_csv(TAB / f"core_collection_{panel.lower()}.csv", index=False)

    pd.DataFrame([dict(panel=panel, n_total=n, recommended_core_size=rec_size,
                        pct_of_panel=round(100 * rec_size / n, 1),
                        pct_diversity_retained=rec["pct_diversity_retained"],
                        full_panel_mean_dist=round(full_mean_dist, 4))]).to_csv(
        TAB / f"core_collection_summary_{panel.lower()}.csv", index=False)

    underused = core_df[~core_df["in_core"]].sort_values("core_rank")
    underused.to_csv(TAB / f"underused_accessions_{panel.lower()}.csv", index=False)
    print(f"  {len(underused)} accessions outside the core -> underused_accessions_{panel.lower()}.csv "
          f"(pre-breeding priority list, ordered by how early they'd re-enter a larger core)")

print("\n-> tables/core_collection_{set}.csv + core_collection_curve_{set}.csv + "
      "core_collection_summary_{set}.csv + underused_accessions_{set}.csv")
