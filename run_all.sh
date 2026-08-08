#!/usr/bin/env bash
# One-command reproduction of the population-structure / diversity / genomic-vs-phenomic
# concordance analysis for the two disjoint rice panels (Set 1 WGS, Set 2 50K array).
# Usage:  bash run_all.sh          (edit the CONFIG block first)
set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"; cd "$ROOT"; export COMPENDIUM_ROOT="$ROOT"

# ============================= CONFIG (edit as needed) =============================
PY="../.venv/Scripts/python.exe"                       # venv python
RS="C:/Program Files/R/R-4.4.3/bin/Rscript.exe"        # R 4.4.3
# Stage toggles (true/false)
RUN_PREP=true      # [0] QC + LD-prune Set1; HapMap->PLINK Set2 -> analysis/data/interim
RUN_PCA=true       # [1] PCA per panel (SNPRelate)
RUN_ADMIX=true     # [2] admixture K-scan (LEA::snmf; ADMIXTURE-equivalent, Windows-native)
RUN_TREE=true      # [3] IBS/kinship + NJ tree
RUN_FST=true       # [4] pairwise Fst + AMOVA by assigned subpopulation
RUN_DIV=true       # [5] diversity: pi, He, Ho, PIC, MAF spectrum, inbreeding F, Ne
RUN_LD=true        # [6] LD decay (r^2 vs distance) -> half-decay distance
RUN_ANCHOR=true    # [7] assign 3K-RGP subpopulations (XI/GJ/cA/cB)
RUN_CONC=true      # [8] PILLAR B: genomic<->phenomic<->trait concordance (Mantel/Procrustes/RV/ARI)
RUN_CORE=true      # [8b] PILLAR C: core/mini-core collection (Python M-strategy; corehunter/rJava crashes)
RUN_FIGS=true      # [9] figures
RUN_TABLES=true    # [10] summary tables
RUN_AUDIT=true     # [11] QC gate
# ---- flagship-upgrade stages (2026-08-08) ----
RUN_ROH=true       # [13] runs of homozygosity + F_ROH
RUN_NE=true        # [14] LD-based Ne (Sved/Corbin, selfing-adjusted)
RUN_WCFST=true     # [15] Weir-Cockerham theta + block-bootstrap CI
RUN_HWE=true       # [16] HWE / Fis / selfing rate
RUN_HAPLO=true     # [17] Beagle phasing + rehh iHS/XP-EHH (Set1; LONG)
RUN_OUTLIER=true   # [18] pcadapt outlier scan
RUN_SFS=true       # [19] folded SFS + Stairway Plot 2 demography (Set1; LONG)
RUN_F3=true        # [20] Patterson f3 admixture tests
RUN_BOOT=true      # [21] NJ bootstrap + DAPC (PowerShell/native shell ONLY: SVD rule)
RUN_RICH=true      # [22] rarefied allelic richness + private alleles
RUN_DIVFIX=true    # [23] diversity recompute on missingness-only SNPs (no MAF filter)
RUN_ROBUST=true    # [24] Pillar B robustness (classifier null, MMRR, sensitivity, FDR)
RUN_LDCHK=true     # [25] LD threshold conventions + thinning-bias check
RUN_AMOVA2=true    # [26] two-level AMOVA (3K-RGP hierarchy)
RUN_FIGS2=true     # [28/29/35] final-scheme figures (Fig06/07 + SuppFig03-05; Fig01; composites Fig02/04/05/08)
RUN_SRCDATA=true   # source data workbook
# ==================================================================================

S=analysis/scripts
[ "$RUN_PREP"   = true ] && { echo "[0] Prep genotypes (QC + LD-prune Set1; HapMap->PLINK Set2)"; "$PY" $S/01_prep_genotypes.py; }
[ "$RUN_PCA"    = true ] && { echo "[1] PCA per panel";                        "$RS" $S/02_pca_structure.R; }
[ "$RUN_ADMIX"  = true ] && { echo "[2] Admixture K-scan (LEA snmf, elbow+global-min K)"; "$RS" $S/03_admixture.R; echo "[2b] UMAP embedding (3rd triangulation method)"; "$PY" $S/03b_umap_structure.py; echo "[2c] Structure consensus (ARI/NMI: PCA vs admixture vs UMAP)"; "$PY" $S/03c_structure_consensus.py; }
[ "$RUN_TREE"   = true ] && { echo "[3] Kinship (GRM) + IBS distance";       "$RS" $S/04_kinship_tree.R; echo "[3b] NJ tree (separate process: SNPRelate+ape crash together)"; "$RS" $S/04b_nj_tree.R; }
[ "$RUN_FST"    = true ] && { echo "[4] Fst (Nei Gst, Python)";                "$PY" $S/05_fst_amova.py; echo "[4b] AMOVA (pegas)"; "$RS" $S/05b_amova.R; }
[ "$RUN_DIV"    = true ] && { echo "[5] Diversity (pi/He/Ho/PIC/F/Ne)";        "$PY" $S/06_diversity.py; }
[ "$RUN_LD"     = true ] && { echo "[6] LD decay";                            "$PY" $S/07_ld_decay.py; }
[ "$RUN_ANCHOR" = true ] && { echo "[7] Set1 3K-RGP subpopulation labels (authoritative, parsed)"; "$PY" $S/08a_parse_set1_subpop.py; echo "[7b] Global 3K-RGP anchor PCA (Set2 not anchored yet, see script header)"; "$RS" $S/08b_anchor_pca.R; }
[ "$RUN_CONC"   = true ] && { echo "[8] Pillar B: genomic<->phenomic<->trait concordance (Python: Mantel/partial-Mantel/Procrustes/RF/feature-attribution/confounder)"; "$PY" $S/09_phenomic_concordance.py; }
[ "$RUN_CORE"   = true ] && { echo "[8b] Pillar C: core/mini-core collection (Python M-strategy max-min diversity)"; "$PY" $S/09b_core_collection.py; }
[ "$RUN_FIGS"   = true ] && { echo "[9] Figures";                              "$PY" $S/10_figures.py; }

# ---- flagship-upgrade stages (13-26; run before the table/figure/audit collectors) ----
[ "$RUN_ROH"    = true ] && { echo "[13] ROH / F_ROH";                        "$PY" $S/13_roh.py; }
[ "$RUN_HWE"    = true ] && { echo "[16] HWE / selfing rate";                 "$PY" $S/16_hwe_selfing.py; }
[ "$RUN_NE"     = true ] && { echo "[14] LD-based Ne (needs 16)";             "$PY" $S/14_ne_ldbased.py; }
[ "$RUN_WCFST"  = true ] && { echo "[15] Weir-Cockerham theta + CI";          "$PY" $S/15_fst_wc.py; }
[ "$RUN_DIVFIX" = true ] && { echo "[23] Diversity recompute (no MAF filter)"; "$PY" $S/23_diversity_recompute.py; }
[ "$RUN_OUTLIER" = true ] && { echo "[18] pcadapt outlier scan";              "$RS" $S/18_selection_outliers.R; }
[ "$RUN_HAPLO"  = true ] && { echo "[17] Beagle phasing + rehh scans (LONG)"; "$RS" $S/17_phase_haplo_scans.R; }
[ "$RUN_SFS"    = true ] && { echo "[19] SFS + Stairway demography (LONG)";   "$PY" $S/19_sfs_demography.py; }
[ "$RUN_F3"     = true ] && { echo "[20] Patterson f3 tests";                 "$PY" $S/20_fstats.py; }
[ "$RUN_BOOT"   = true ] && { echo "[21] NJ bootstrap + DAPC";                "$RS" $S/21_nj_boot_dapc.R; }
[ "$RUN_RICH"   = true ] && { echo "[22] Allelic richness + private alleles"; "$PY" $S/22_richness.py; }
[ "$RUN_ROBUST" = true ] && { echo "[24] Pillar B robustness";                "$PY" $S/24_robustness.py; }
[ "$RUN_LDCHK"  = true ] && { echo "[25] LD thresholds + thinning check";     "$PY" $S/25_ld_check.py; }
[ "$RUN_AMOVA2" = true ] && { echo "[26] Two-level AMOVA";                    "$RS" $S/26_amova2.R; }
[ "$RUN_FIGS2"  = true ] && { echo "[28] Flagship figures";                   "$PY" $S/28_figures_flagship.py; echo "[29] Fig01 design schematic (house standard)"; "$PY" $S/29_fig01_design.py; echo "[36] Mechanism/replication source tables (Fig09/Fig10)"; "$PY" $S/36_mechanism_visuals.py; echo "[35] Composite high-density figures"; "$PY" $S/35_figures_composite.py; }
# ---- review-response stages (2026-08-08, F1-F8 fixes) ----
RUN_ROBUST=true
[ "$RUN_ROBUST" = true ] && { echo "[30] Confound + label-robustness battery (LONG)"; "$PY" $S/30_confound_robustness.py; echo "[31] Demography repair (F5)"; "$PY" $S/31_demography_repair.py; echo "[32] Scan-overlap rotation null (F6)"; "$PY" $S/32_overlap_null.py; echo "[33] Reference validation + Procrustes reconciliation (F7)"; "$PY" $S/33_validate_reference.py; echo "[37] Size-leakage + classifier uncertainty (LONG)"; "$PY" $S/37_mechanism_stature_tests.py; echo "[38] Mechanism forensics: K-forensics, residualised clf, temporal Mantel, Pst-Fst, FDR ledger"; "$PY" $S/38_mechanism_forensics.py; echo "[38b] Residualised-classifier permutation null (LONG)"; "$PY" $S/38b_residualised_perm.py; }

[ "$RUN_TABLES" = true ] && { echo "[10] Tables (1-7, S1-S16)";               "$PY" $S/11_tables.py; }
[ "$RUN_SRCDATA" = true ] && { echo "[SD] Source Data workbook";              "$PY" $S/build_source_data.py; }
[ "$RUN_AUDIT"  = true ] && { echo "[11] QC gate";                            "$PY" $S/12_audit.py; }

echo "DONE. Outputs: analysis/data/interim/ (PLINK/GDS), analysis/results/{figures,tables,supp_tables,source_data}/"