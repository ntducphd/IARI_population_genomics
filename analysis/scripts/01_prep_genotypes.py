#!/usr/bin/env python
# 01_prep_genotypes.py — QC + LD-prune both disjoint panels into clean PLINK bfiles for structure analysis.
#   Set1 (WGS, PLINK bfile) : MAF>0.05, geno<0.10, then LD-prune (indep-pairwise 50 5 0.2)
#   Set2 (50K, HapMap)      : HapMap -> VCF (count the ALT allele) -> PLINK, same QC + prune
# Outputs (analysis/data/interim/): {Set1,Set2}_qc.{bed,bim,fam} (full, for diversity/LD) and
#   {Set1,Set2}_pruned.{bed,bim,fam} (unlinked, for PCA/admixture) + prep_summary.csv
import subprocess, csv
from paths import PLINK, SET1_BED_1M, SET2_HMP, INTERIM, TAB

CHRSET = ["--chr-set", "12", "no-xy"]
def plink(*args):
    cmd = [str(PLINK), *[str(a) for a in args], *CHRSET, "--allow-extra-chr", "--silent"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("PLINK ERROR:\n", r.stderr[-1500:]); raise SystemExit(1)
    return r

def counts(stem):
    import pathlib
    b = pathlib.Path(f"{stem}.bim"); f = pathlib.Path(f"{stem}.fam")
    return (sum(1 for _ in open(b)) if b.exists() else 0,
            sum(1 for _ in open(f)) if f.exists() else 0)

def qc_and_prune(in_stem, tag, vcf=False):
    qc = INTERIM / f"{tag}_qc"; pr = INTERIM / f"{tag}_pruned"
    src = ["--vcf", in_stem] if vcf else ["--bfile", in_stem]
    plink(*src, "--maf", 0.05, "--geno", 0.10, "--make-bed", "--out", qc)
    plink("--bfile", qc, "--indep-pairwise", 50, 5, 0.2, "--out", INTERIM / f"{tag}_ldp")
    plink("--bfile", qc, "--extract", INTERIM / f"{tag}_ldp.prune.in", "--make-bed", "--out", pr)
    return qc, pr

# ---- Set2 HapMap -> VCF (ALT-allele dosage) ----
def hmp_to_vcf(hmp, out_vcf):
    with open(hmp) as fh:
        header = fh.readline().rstrip("\n").split(",")
        samples = header[11:]
        MISS = {"NN", "--", "NA", "", "??"}
        with open(out_vcf, "w") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
            n = 0
            for line in fh:
                p = line.rstrip("\n").split(",")
                snp, alleles, chrom, pos = p[0], p[1], p[2], p[3]
                if "/" not in alleles:  # skip malformed
                    continue
                ref, alt = alleles.split("/")[0], alleles.split("/")[1]
                gts = []
                for g in p[11:]:
                    g = g.strip().upper()
                    if g in MISS or len(g) != 2:
                        gts.append("./."); continue
                    a1, a2 = g[0], g[1]
                    def code(a): return "0" if a == ref else ("1" if a == alt else ".")
                    c1, c2 = code(a1), code(a2)
                    gts.append(f"{c1}/{c2}" if "." not in (c1, c2) else "./.")
                out.write(f"{chrom}\t{pos}\t{snp}\t{ref}\t{alt}\t.\t.\t.\tGT\t" + "\t".join(gts) + "\n")
                n += 1
    return n, len(samples)

INTERIM.mkdir(parents=True, exist_ok=True)
rows = []

# Set1
print("[Set1] QC + LD-prune (WGS bfile)")
raw1 = counts(str(SET1_BED_1M))
qc1, pr1 = qc_and_prune(str(SET1_BED_1M), "Set1")
rows.append(("Set1", "WGS 1.01M", raw1[0], counts(str(qc1))[0], counts(str(pr1))[0], counts(str(pr1))[1]))

# Set2
print("[Set2] HapMap -> VCF -> QC + LD-prune (50K array)")
vcf2 = INTERIM / "Set2.vcf"
nsnp2, nsamp2 = hmp_to_vcf(SET2_HMP, vcf2)
print(f"       hmp->vcf: {nsnp2} SNPs x {nsamp2} accessions")
qc2, pr2 = qc_and_prune(str(vcf2), "Set2", vcf=True)
rows.append(("Set2", "50K array", nsnp2, counts(str(qc2))[0], counts(str(pr2))[0], counts(str(pr2))[1]))

# summary
out_csv = TAB / "prep_summary.csv"
with open(out_csv, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["panel", "platform", "snps_raw", "snps_qc", "snps_pruned", "n_accessions"])
    w.writerows(rows)
print("\n=== PREP SUMMARY ===")
print(f"{'panel':6}{'platform':14}{'raw':>10}{'qc':>10}{'pruned':>10}{'n':>6}")
for r in rows:
    print(f"{r[0]:6}{r[1]:14}{r[2]:>10}{r[3]:>10}{r[4]:>10}{r[5]:>6}")
print(f"-> {out_csv}")
