#!/usr/bin/env python
"""19_sfs_demography.py — [stage 19] folded site-frequency spectrum + Stairway Plot 2 demographic
history, Set 1 only.

Motivation: the genome-history layer. Set 1 is
WGS-derived so its SFS is usable; Set 2 (50K chip) is EXCLUDED — array SNPs are ascertained on a
discovery panel, which distorts the SFS by design (stated in the manuscript as a teaching point,
not a footnote).

Method decisions (each stated in Methods):
  * Input = the missingness-only-filtered Set 1 genotypes (stage 23) — an SFS computed after a
    MAF filter would be missing its low-frequency classes entirely.
  * HAPLOIDISATION: with near-complete selfing (F ~ 0.93, stage 16) the two alleles of one
    accession are pseudo-replicates; standard practice for inbred-line panels is one allele per
    accession per site (random draw, seed 42) -> n = 150 haploid samples.
  * FOLDED spectrum (no reliable ancestral-allele polarisation without outgroup alignment).
  * Stairway Plot 2 (Liu & Fu 2020) with the blueprint's four random break points; mutation rate
    and generation time are SCALING constants only: mu = 7.0e-9 /bp/generation (order of the
    values used in rice demography literature; the manuscript reports the trajectory in scaled
    and absolute time with this stated) and 1 generation/year. L = callable length proxied by the
    reference genome span covered by the SNP set (sum of per-chromosome max positions) — an
    approximation; stated.
  * Bootstrap inputs reduced to NINPUT=30 (compute budget); CIs reported as exploratory.

Outputs:
  tables/sfs_folded_set1.csv        — allele count class k (1..75), n_snps
  tables/stairway_ne_set1.csv       — year, Ne_median, Ne_2.5%, Ne_97.5%  (parsed final summary)
  interim/stairway/                  — Stairway working directory (regenerable)
"""
import subprocess
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import allel

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import INTERIM, TAB, ROOT

TOOLS = ROOT / "analysis/tools"
SPDIR = TOOLS / "stairway_plot_v2.1.2"
WORK = INTERIM / "stairway"
MU = 7.0e-9
GEN_YEARS = 1
NINPUT = 30
RNG = np.random.default_rng(42)


def folded_sfs():
    cs = allel.read_vcf(str(INTERIM / "Set1_nomaf.vcf"),
                        fields=["variants/CHROM", "variants/POS", "calldata/GT"])
    gt = allel.GenotypeArray(cs["calldata/GT"])       # (loci, 150, 2)
    chrom, pos = cs["variants/CHROM"], cs["variants/POS"]
    # haploidise: one random allele per accession per site; missing stays missing
    pick = RNG.integers(0, 2, size=gt.shape[:2])
    hap = np.take_along_axis(np.asarray(gt), pick[..., None], axis=2)[..., 0]
    hap = np.ma.masked_equal(hap, -1)
    n = hap.shape[1]
    alt = hap.sum(axis=1).filled(-1)
    called = (~hap.mask).sum(axis=1)
    # keep fully-called sites so every SNP contributes to the same n=150 spectrum
    keep = (called == n) & (alt >= 0)
    alt = alt[keep]
    minor = np.minimum(alt, n - alt)
    sfs = np.bincount(minor, minlength=n // 2 + 1)[1:]           # classes 1..n/2
    L = int(pd.DataFrame({"c": chrom, "p": pos}).groupby("c")["p"].max().sum())
    return sfs, n, L, int(keep.sum())


def write_blueprint(sfs, n, L):
    WORK.mkdir(parents=True, exist_ok=True)
    nseq = n
    bp = WORK / "set1.blueprint"
    lines = [
        "popid: Set1",
        f"nseq: {nseq}",
        f"L: {L}",
        "whether_folded: true",
        "SFS: " + " ".join(str(int(x)) for x in sfs),
        "smallest_size_of_SFS_bin_used_for_estimation: 1",
        f"largest_size_of_SFS_bin_used_for_estimation: {nseq // 2}",
        "pct_training: 0.67",
        # four random break points per Stairway convention: (n-2)/4, (n-2)/2, 3(n-2)/4, n-2
        f"nrand: {(nseq - 2) // 4} {(nseq - 2) // 2} {3 * (nseq - 2) // 4} {nseq - 2}",
        # Stairbuilder's path handling chokes on absolute Windows paths -> keep everything
        # relative to the working directory (stairway_plot_es is copied in by main()).
        "project_dir: set1_out",
        "stairway_plot_dir: stairway_plot_es",
        f"ninput: {NINPUT}",
        "random_seed: 42",
        f"mu: {MU}",
        f"year_per_generation: {GEN_YEARS}",
        "plot_title: Set1",
        "xrange: 0,0",
        "yrange: 0,0",
        "xspacing: 2",
        "yspacing: 2",
        "fontsize: 12",
    ]
    bp.write_text("\n".join(lines) + "\n")
    return bp


def main():
    if not SPDIR.exists():
        with zipfile.ZipFile(TOOLS / "stairway_plot_v2.zip") as z:
            z.extractall(TOOLS)
    WORK.mkdir(parents=True, exist_ok=True)
    if not (WORK / "stairway_plot_es").exists():
        import shutil
        shutil.copytree(SPDIR / "stairway_plot_es", WORK / "stairway_plot_es")

    sfs, n, L, n_snps = folded_sfs()
    pd.DataFrame({"minor_allele_count": np.arange(1, len(sfs) + 1), "n_snps": sfs}).to_csv(
        TAB / "sfs_folded_set1.csv", index=False)
    print(f"[19] folded SFS: {n_snps} fully-called SNPs, n={n} haploids, L={L}")
    print(f"[19] singleton share: {sfs[0] / sfs.sum():.3f}")

    bp = write_blueprint(sfs, n, L)
    r = subprocess.run(["java", "-cp", str(SPDIR / "stairway_plot_es"), "Stairbuilder",
                        str(bp)], capture_output=True, text=True, cwd=WORK)
    if r.returncode != 0:
        print("Stairbuilder failed:", r.stderr[-2000:])
        raise SystemExit(1)
    # Stairbuilder writes "<blueprint>.sh" beside the blueprint
    print(f"[19] running Stairway batch ({NINPUT} inputs)...")
    # Stairbuilder emits a .bat on Windows and a .sh elsewhere -> run whichever exists
    bat, sh = Path(str(bp) + ".bat"), Path(str(bp) + ".sh")
    if bat.exists():
        r2 = subprocess.run(["cmd", "/c", str(bat)], capture_output=True, text=True, cwd=WORK)
    else:
        git_bash = Path("C:/Program Files/Git/bin/bash.exe")
        bash = str(git_bash) if git_bash.exists() else "bash"
        r2 = subprocess.run([bash, sh.name], capture_output=True, text=True, cwd=WORK)
    if r2.returncode != 0:
        print("Stairway batch failed:", (r2.stderr or r2.stdout)[-2000:])
        raise SystemExit(1)

    final = list((WORK / "set1_out").glob("*final.summary"))
    if not final:
        print("[19] no final.summary found"); raise SystemExit(1)
    fs = pd.read_csv(final[0], sep=r"\s+")
    out = fs.rename(columns={"year": "year", "Ne_median": "ne_median",
                             "Ne_2.5%": "ne_lo95", "Ne_97.5%": "ne_hi95"})
    keep = [c for c in ["mutation_per_site", "n_estimation", "year", "ne_median",
                        "ne_lo95", "ne_hi95"] if c in out.columns]
    out[keep].to_csv(TAB / "stairway_ne_set1.csv", index=False)
    print(f"[19] done -> stairway_ne_set1.csv ({len(out)} steps); "
          f"most-recent Ne_median = {out['ne_median'].iloc[0]:.0f}")


if __name__ == "__main__":
    main()
