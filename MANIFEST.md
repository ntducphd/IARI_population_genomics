# MANIFEST — reproducibility map (float → script → inputs → outputs)

Reproducibility contract for this compendium: for every planned figure/table it names the
producing script and the data it consumes. With `run_all.sh` (the master script) a reader regenerates
every result in manuscript order. The *analysis is the workflow*. Scripts marked **[to build]** are the
committed plan (built in numbered order); `paths.py` / `_paths.R` are the single source of truth for paths.

## Reproduce in one command
```bash
bash run_all.sh          # stages [0]–[11]; edit the CONFIG block to toggle stages
```
Outputs land in `analysis/data/interim/` (PLINK/GDS) and `analysis/results/{figures,tables,supp_tables}/`.

## Pipeline execution order (`run_all.sh`)
| Stage | Script | Role |
|-------|--------|------|
| [0] | `01_prep_genotypes.py` **[built, ran]** | Set1: PLINK QC (MAF>0.05, geno<0.1) + LD-prune; Set2: HapMap→PLINK; harmonise IDs → `data/interim/{Set1,Set2}_qc[.prune].{bed,bim,fam}` + `.gds`. Real result: `tables/prep_summary.csv` (Set1 502,675 QC/24,370 pruned/150 acc; Set2 31,565 QC/2,942 pruned/147 acc). |
| [1] | `02_pca_structure.R` **[built, ran]** | SNPRelate PCA per panel → `tables/pca_{set}.csv` (PCs + %var), eigenvalues. Real result: `tables/pca_variance.csv` (dominant k-means cluster 50% Set1, 53.1% Set2 -- the narrow-elite-base premise, confirmed). |
| [2] | `03_admixture.R` **[built]** | Set1/Set2: PLINK->VCF->LEA::vcf2geno->snmf K=1..12, 5 reps; K selected by BOTH elbow (headline) and global-min (transparency) cross-entropy → `tables/admixture_{set}_Q.csv`, `admixture_{set}_cv.csv`, `admixture_bestK.csv` |
| [2b] | `03b_umap_structure.py` **[built]** | 3rd triangulation method (DAPC dropped -- adegenet::read.PLINK crashes R on this machine, see the environment-audit record): PLINK --recode A -> mean-impute -> UMAP(2D) → `tables/umap_{set}.csv` |
| [2c] | `03c_structure_consensus.py` **[built]** | ARI/NMI concordance between PCA k-means, sNMF admixture (argmax Q), and UMAP k-means, at both k=3 and the elbow K → `tables/structure_consensus_{set}.csv`, `structure_consensus_summary.csv`. Same ARI/NMI toolkit that Pillar B (stage 8) reuses for genomic<->phenomic concordance. |
| [3] | `04_kinship_tree.R` **[built, ran]** | GRM kinship (GCTA method) + IBS distance (SNPRelate+gdsfmt only) → `tables/kinship_{set}.csv`, `ibs_dist_{set}.csv`, `interim/ibs_dist_{set}.rds`. Real: Set1 mean IBS-dist 0.310 (range 0.008-0.396), Set2 mean 0.309 (range 0.004-0.457). |
| [3b] | `04b_nj_tree.R` **[built, ran]** | NJ tree from the IBS distance (ape only, separate process -- SNPRelate+gdsfmt+ape crash together on this machine, see script header) → `interim/nj_{set}.nwk`, `tables/nj_{set}_tips.csv`. Bootstrap deferred (STRETCH: IQ-TREE). |
| [4] | `05_fst_amova.py` **[built, ran]** | Pairwise + global Fst between admixture clusters (Nei's Gst Hs/Ht, computed directly in Python -- hierfstat abandoned, see script header) + permutation test (999 perms) on the LD-pruned set → `tables/fst_pairwise_{set}.csv`, `fst_global_{set}.csv`. Real: Set1 global Fst=0.330 (p=0.001), Set2 global Fst=0.294 (p=0.001), pairwise range 0.05-0.39, comparable to Zhao et al. 2011's 0.23-0.53 across 5 rice subpopulations. |
| [4b] | `05b_amova.R` **[built, ran]** | Hierarchical AMOVA (pegas only, on the same IBS distance + admixture-cluster grouping) → `tables/amova_{set}.csv`. Real: 56.4% (Set1) / 60.9% (Set2) of variance among admixture clusters -- consistent with the Fst estimates. |
| [5] | `06_diversity.py` **[built, ran]** | π, θ_W (windowed, 1Mb, per-chromosome then site-count-weighted genome-wide combine -- naive cross-chromosome concatenation crashes allel.sequence_diversity, fixed), He/Ho/PIC/MAF spectrum (PLINK --freq/--het), Tajima's D. Ne (LD-based) explicitly NOT computed (deferred, would need a correctly-reimplemented estimator). Real results: F=0.926 (Set1) / 0.787 (Set2) -- very high inbreeding, expected for a self-pollinating rice breeding panel; Tajima's D +3.15/+2.76 -- consistent with the narrow-elite-base/breeding-bottleneck story → `tables/diversity_summary.csv`, `maf_spectrum_{set}.csv`, `diversity_windows_{set}.csv` |
| [6] | `07_ld_decay.py` **[built, ran]** | PLINK `--r2` (Set1 thinned 10% first -- untinned would emit ~10^8 pairs, intractable; Set2 used directly, already sparse) → bin r² vs distance (25kb bins) → half-decay distance → `tables/ld_decay_{set}.csv`, `ld_decay_summary.csv`. Real: Set1 half-decay ~562 kb (r2@12.5kb=0.213, thinning likely inflates this vs full density), Set2 ~212 kb (r2@12.5kb=0.425, close to the ~196 kb DRC-rice DArTseq literature reference). |
| [7] | `08a_parse_set1_subpop.py` **[built, ran]** | Parses the AUTHORITATIVE 3K-RGP/SNP-Seek subpopulation labels for Set1 already on disk (`data/raw/genotype/.../mylists-634211649175051595.txt`) rather than re-deriving them → `tables/subpop_assignment_set1.csv`, `subpop_capture_set1.csv`. Real: aus 34.0% + indx 24.0% + ind2 18.7% = 76.7% of Set1 in just 3 indica-related groups -- strong narrow-elite-base evidence from ground-truth labels, not our own clustering. |
| [7b] | `08b_anchor_pca.R` **[built, ran]** | PCA of Set1 projected into the merged 3024-accession 3K-RGP+Set1 reference panel (`pruned_v2.1`, thinned ~10% first -- untinned BED->GDS did not finish in 19 min, killed) → `tables/anchor_pca_global.csv`, `anchor_pca_variance.csv`. Real, textbook-consistent result: indica-related groups (ind1A/B, ind2/3, indx) all PC1>0, japonica-related groups (japx/subtrop/temp/trop) all PC1<0, aus separates on PC2 -- matches the classic indica/japonica split (Wang 2018, Huang 2012). **Set2 (50K array) NOT anchored**: no pre-merged 3K-RGP+Set2 genotype file found on disk, only a common-SNP-position list without the actual reference genotypes -- documented gap, not fabricated; would need additional data preparation to close. |
| [8] | `09_phenomic_concordance.py` **[built, ran]** | PILLAR B (the paper's central analysis): Mantel/partial-Mantel + Procrustes/PROTEST + RF classification (phenomic->admixture cluster) + feature-type attribution (Colour/NIR/Size) + structure-as-confounder test, ALL implemented directly in Python (not R -- see script header) → `tables/concordance_{mantel,procrustes,classification,feature_attribution,confounder}_{set}.csv`. Real, cross-panel-replicated result: genomic<->phenomic Mantel r=0.160 (Set1, p=0.001) / r=0.087 (Set2, p=0.004) -- imaging carries a significant, independently-replicated genetic-structure signal; Procrustes p=0.001 both panels; RF phenomic->cluster beats majority baseline 2x (Set1: 45.2% vs 22.6%; Set2: 29.9% vs 19.1%); Size/morphology features carry the strongest genomic signal in both panels, Colour is weakest (not significant in Set2, p=0.348); genomic structure is only a MODEST confounder of the phenomic-trait relationship (+2.5-5.1% change when partialled out) -- reassuring for the companion phenomic-selection study's claims. Uses the case-insensitive genotype-phenotype bridge rule from data/README.md (naive exact-match recovers only 109/147 for Set2, the documented known-incomplete number; fixed, now 147/147). |
| [8b] | `09b_core_collection.py` **[built, ran]** | PILLAR C breeding resource: greedy max-min-diversity ("M strategy", Schoen & Brown 1993) core-collection selection in Python -- corehunter (R/rJava) installs but crashes at `.jinit()`, reproducibly (documented environment audit). Recommends the conventional ~10%-of-panel core size (Frankel & Brown 1984), NOT a naive diversity-threshold search (max-min sampling inflates retained-diversity above 100% at very small sizes by construction, which degenerates a threshold rule to an unusably tiny core -- caught and fixed) → `tables/core_collection_{set}.csv`, `core_collection_curve_{set}.csv`, `core_collection_summary_{set}.csv`, `underused_accessions_{set}.csv`. Real: 15-accession (~10%) core retains 101.5% (Set1) / 107.7% (Set2) of full-panel mean pairwise diversity. |
| [9] | `10_figures.py` **[built, ran; re-scoped 2026-08-08 to the 8-figure composite scheme]** | Under the final scheme this stage produces **Fig03_tree_kinship** (NJ trees + GRM heatmaps, tips coloured by admixture cluster) plus SuppFig01 (UMAP+ARI consensus) and SuppFig02 (3K-RGP anchor); its former single-topic mains (PCA, admixture, Fst+diversity, LD, concordance, core-curve supp) now render as `_archive_*` files (kept for teaching, outside the submission bundle glob) — their content lives on inside the composite Figs 2/4/5/8 built by stage [35]. Shared Okabe-Ito cluster palette + C0.. labels via `figstyle.CLUSTER_COL`. |
| [10] | `11_tables.py` **[built, ran]** | Table 1 (panel/K-selection/diversity), Table 2 (Fst+AMOVA), Table 3 (LD decay), Table 4a-e (Pillar B: Mantel/Procrustes/classification/feature-attribution/confounder), Table 5 (core collection) → `results/tables/Table_*.{md,csv}`. Supp Tables S1-S7 (prep detail, subpop capture, full pairwise Fst, MAF spectrum, structure consensus, core-collection curve, underused accessions) → `results/supp_tables/`. |
| [11] | `12_audit.py` **[built, ran]** | QC gate: figure/table existence (PNG+PDF pairs), p-value/correlation range sanity, headline-result cross-panel replication check, NaN scan. **PASSED: 34 checks, 1 known warning (Set2 3K-RGP gap, documented), 0 errors.** |

## Float → producing script → inputs → output (final 10-figure scheme; supp figures S1–S5; all paths under `analysis/`)
| Float | Script | Key inputs | Output |
|-------|--------|-----------|--------|
| Fig 1 · study design + workflow schematic | `29_fig01_design.py` | panel metadata | `results/figures/main/Fig01_design.*` |
| Fig 2 · structure & ancestry (6p: PCA a-b, cross-entropy c-d, Q-bars e-f) | `35_figures_composite.py` ← `02,03` | `tables/pca_set{1,2}.csv`, `admixture_*_{Q,cv}.csv` | `Fig02_structure.*` |
| Fig 3 · NJ trees + GRM kinship (4p) | `10_figures.py` ← `04` | `interim/nj_*.nwk`, `tables/kinship_*.csv` | `Fig03_tree_kinship.*` |
| Fig 4 · differentiation + diversity + LD (4p: Gst a-b, He/Ho/F c, LD decay d) | `35_figures_composite.py` ← `05,06,07,23` | `tables/fst_*.csv`, `diversity_*.csv`, `ld_decay_*.csv` | `Fig04_diff_diversity_ld.*` |
| Fig 5 · concordance + mechanism (6p: Mantel CI, Procrustes band, classifier nulls, feature families, stature collapse, subset classifiers) | `35_figures_composite.py` ← `09,24,30,33` | `tables/concordance_*`, `confound_mantel_*`, `classifier_baselines_*`, `mantel_bootstrap_*`, `robustness_procrustes_*` | `Fig05_concordance_mechanism.*` |
| Fig 6 · genome history (4p: ROH classes, F_ROH by cluster, LD-Ne, SFS) | `28_figures_flagship.py` ← `13,14,19` | `tables/roh_*`, `ne_trajectory_*`, `sfs_folded_set1` | `Fig06_genome_history.*` |
| Fig 7 · selection scans (4p: pcadapt ×2, iHS, XP-EHH) | `28_figures_flagship.py` ← `17,18` | `tables/pcadapt_*`, `ihs_set1`, `xpehh_set1` | `Fig07_selection.*` |
| Fig 8 · application (3p: core curves, private alleles, rarefied richness) | `35_figures_composite.py` ← `08b,22` | `tables/core_collection_curve_*`, `richness_*` | `Fig08_application.*` |
| Fig 9 · mechanism made visible (6p: height-by-cluster a-b, distance-pair hexbins c-d, external confusion e, robustness ladder f) | `35_figures_composite.py` ← `36,30` | `tables/height_by_cluster_*`, `distance_pairs_*`, `external_confusion_set1`, `classifier_baselines_*` | `Fig09_mechanism_visible.*` |
| Fig 10 · two platforms, one biology (3p: replication scoreboard a, MAF spectra b, SNP density c) | `35_figures_composite.py` ← `36,32` | `tables/replication_scoreboard`, `maf_spectrum_*`, `snp_density_by_chrom` | `Fig10_platforms.*` |
| Supp S1 UMAP+consensus, S2 3K-RGP anchor | `10_figures.py` | — | `figures/supp/SuppFig0{1,2}_*` |
| Supp S3 Stairway, S4 f3, S5 DAPC | `28_figures_flagship.py` ← `19,20,21` | — | `figures/supp/SuppFig0{3,4,5}_*` |
| Table 1 · panel/diversity summary | `11_tables.py` ← `06` | `tables/diversity_*.csv` | `tables/Table_1_diversity.{md,csv}` |
| Table 2 · F_ST / AMOVA | `11_tables.py` ← `05` | `tables/fst_*.csv`, `amova_*.csv` | `tables/Table_2_fst_amova.{md,csv}` |

## Cross-paper dependencies (fixed relative paths; read-only)
- `01_prep_genotypes.py` ← `../data/raw/genotype/{Subset1_150Geno_1M/Genotypes/150genotypes.*, Subset2_147Geno/Genotypes/147SNPgenoypes.hmp.csv}`, `../manuscript_6_gwas_nue/scripts/Plink1.9/plink.exe`
- `08b_anchor_pca.R` ← `../data/raw/genotype/Subset1_150Geno_5.2M/Wanget al.2018/3K-HDRA-snp-comm-miss5pc.txt` (+ 3K-RGP reference genotypes)
- `09_phenomic_concordance.py` ← `../manuscript_7_phenomic_selection/analysis/data/input/cohort_set{1,2}.csv` (204 phenomic features + traditional traits, ID-aligned)

## Flagship-upgrade stages (13-26, added 2026-08-08)

| Stage | Script | Role | Key outputs (analysis/results/tables/) |
|---|---|---|---|
| [13] | 13_roh.py | ROH + F_ROH (PLINK --homozyg, per-density params) | roh_summary, roh_indiv_*, roh_length_classes_*, roh_cluster_summary_* |
| [14] | 14_ne_ldbased.py | LD-based Ne (Sved/Weir-Hill/Corbin; selfing c*=c(1-F)) | ne_summary, ne_trajectory_* |
| [15] | 15_fst_wc.py | Weir-Cockerham theta + 1-Mb block-bootstrap CI + cluster sizes | fst_wc_global_*, fst_wc_pairwise_*, cluster_sizes_* |
| [16] | 16_hwe_selfing.py | HWE exact tests, Fis, selfing rate s=2F/(1+F) | hwe_summary, fis_per_accession_* |
| [17] | 17_phase_haplo_scans.R | Beagle 5.5 phasing -> rehh iHS + XP-EHH (Set1) | ihs_set1, xpehh_set1, selection_regions_set1, haplo_scan_summary |
| [18] | 18_selection_outliers.R | pcadapt genome scan (BH q<0.05) | pcadapt_outliers_*, pcadapt_summary |
| [19] | 19_sfs_demography.py | folded SFS (no-MAF, haploidised) + Stairway Plot 2 | sfs_folded_set1, stairway_ne_set1 |
| [20] | 20_fstats.py | Patterson f3 (external 3K-RGP groups, block jackknife) | f3_set1, f3_summary |
| [21] | 21_nj_boot_dapc.R | NJ bootstrap (100x) + DAPC (4th structure method) | nj_boot_support_*, nj_boot_summary, dapc_assign_*, dapc_status |
| [22] | 22_richness.py | rarefied allelic richness + private alleles per cluster | richness_*, richness_summary |
| [23] | 23_diversity_recompute.py | pi/theta_W/Tajima D on missingness-only SNPs (T1-1 fix) | diversity_nomaf_summary, diversity_nomaf_windows_* |
| [24] | 24_robustness.py | classifier permutation null, MMRR, Procrustes/K sensitivity, FDR ledger | robustness_* |
| [25] | 25_ld_check.py | r2=0.2/0.1 threshold decay + empirical thinning check | ld_threshold_summary, ld_thinning_check |
| [26] | 26_amova2.R | two-level AMOVA (macro-group/3K-RGP subpop) | amova2_set1 |
| [28] | 28_figures_flagship.py | Fig06_genome_history, Fig07_selection, SuppFig03-05 (renamed 2026-08-08 from Fig08/09 + SuppFig04-06 in the 8-figure consolidation) | figures/main + figures/supp |
| [SD] | build_source_data.py | Source_Data.xlsx (one sheet per data-bearing main figure) | source_data/Source_Data.xlsx |

New floats (final names): Fig06_genome_history <- [13,14,19]; Fig07_selection <- [17,18]; SuppFig03_stairway <- [19];
SuppFig04_f3 <- [20]; SuppFig05_dapc <- [21]; Table_6_genome_history <- [13,14]; Table_7_selection <- [17,18];
Table 1/2/3 upgrades <- [23,16,13 / 15,26 / 25]; SuppTables S8-S12 <- [20,22,13,24,24].
Third-party tools in analysis/tools/: beagle.jar (Beagle 5.5), stairway_plot_v2.1.2 (fetched by stages; see .gitignore).
| [29] | 29_fig01_design.py | Fig01 study-design schematic built to the house schematic standard (replaces the 2026-08-04 flowchart) | figures/main/Fig01_design.png+pdf |
| [35] | 35_figures_composite.py | High-density composite mains for the final 10-figure scheme: Fig02_structure (6p), Fig04_diff_diversity_ld (4p), Fig05_concordance_mechanism (6p), Fig08_application (3p), Fig09_mechanism_visible (6p), Fig10_platforms (3p). Replaces 34_figures_nature.py, which was DELETED 2026-08-08 together with all _archive_* figure renders (superseded renders removed to avoid stale-name confusion; see git history). | figures/main/Fig{02,04,05,08,09,10}_*.png+pdf |
| [36] | 36_mechanism_visuals.py | Source tables for Fig 9/10 (reuses stage-30 helpers via importlib; external-confusion accuracy ASSERTED equal to the stage-30 headline 58.9% so figure and Table S14 can never drift): height_by_cluster_* + summary (Kruskal-Wallis), distance_pairs_* (10,585/10,731 pairs), external_confusion_set1, replication_scoreboard (12 dimensionless headline metrics gathered from existing summaries — nothing recomputed) | tables/height_by_cluster_*, distance_pairs_*, external_confusion_set1, replication_scoreboard |

## Review-response stages (30-33, added 2026-08-08 after the adversarial review round)

| Stage | Script | Role | Key outputs |
|---|---|---|---|
| [30] | 30_confound_robustness.py | F1/F3/F4: stature-confound Mantel battery, position-proxy tests, height-only/feature-subset classifier baselines, external-3K-RGP-label classifier, K+/-1 targets, soft-Q accuracy, 999-perm nulls, bootstrap CIs | confound_mantel_*, spatial_proxy_*, classifier_baselines_*, external_label_classifier_set1, mantel_bootstrap_* |
| [31] | 31_demography_repair.py | F5: Stairway summary cleaned (Ne columns, underflow rows), within-largest-cluster Tajima D + LD-Ne, Ne interval across F-source x map-density grid | stairway_ne_set1 (clean), within_cluster_history_set1, ne_interval_summary |
| [32] | 32_overlap_null.py | F6: circular-rotation null for two-scan region overlap + exact-SNP overlap; pcadapt saturation + per-chromosome SNP density bookkeeping | overlap_null_summary, pcadapt_saturation, snp_density_by_chrom |
| [33] | 33_validate_reference.py | F7: bespoke Mantel vs vegan::mantel; Procrustes sensitivity recomputed with the stage-09 input pipeline (reconciles the M2 band) | mantel_vegan_validation, robustness_procrustes_* (overwritten) |
| [37] | 37_mechanism_stature_tests.py | Size-leakage CV-R2 (height index predicted from colour+NIR: ridge 0.72/0.87, RF 0.63/0.80 — the "no size information" claim measured and found FALSE) + 20-repeat CV uncertainty for every classifier cell | sizeleakage_summary, classifier_uncertainty |
| [38] | 38_mechanism_forensics.py | K±1 forensics (byte-identical accuracies = partition coincidences, labels differ for >90% of accessions; primary-protocol folds documented), size-residualised colour+NIR classifier (0.260/0.238), per-imaging-day Mantel (stable, all P<=0.014), Pst-Fst c/h2 sweep (all Pst << theta), FDR ledger extended 16->34 tests (load-bearing marginals survive BH) | kbattery_forensics, residualised_classifier, temporal_mantel, pst_fst_sweep, robustness_fdr_mantel (extended) |
| [38b] | 38b_residualised_perm.py | 999-label-permutation null for the size-residualised classifier (its accuracy is now load-bearing) | residualised_classifier (empirical_p/null_mean added) |

Headline outcomes: two-scan "convergence" is chance-level (enrichment 0.95x, P=0.917; exact-SNP 0.4x
DEPLETED) and is retracted in the manuscript; within-cluster Tajima D = +0.94 (vs +1.14 pooled);
Ne reported as intervals (Set1 305-1267 selfing-adj / 57-94 panmictic); Mantel matches vegan to 3e-5;
Procrustes primary values sit inside the reconciled 2-6-PC band. SuppTables S13-S16 carry these.