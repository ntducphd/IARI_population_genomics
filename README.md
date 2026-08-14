# Canopy phenomics carries genetic population structure: a layered population-genomic portrait of two disjoint rice breeding panels

**Authors:** Nguyen Trung Duc^1,2^ and Dhandapani Raju^1^ (corresponding, dandyman2k6@gmail.com)
^1^Nanaji Deshmukh Plant Phenomics Centre (NDPPC), Division of Plant Physiology,
ICAR-Indian Agricultural Research Institute (IARI), New Delhi, India
^2^Vietnam National University of Agriculture, Hanoi, Vietnam

Reproducible analysis compendium for the manuscript above: all code that produces every
figure, table and supplementary table of the paper, together with the produced result
tables and figure files.

## The analyses

Three hundred rice accessions were split into **two genetically disjoint panels**: Set 1
(150 accessions, whole-genome-sequenced) and Set 2 (147 accessions, 50K array). The panels
share no accessions and were genotyped on different platforms. The pipeline builds a layered
population-genomic portrait:

- **Structure**: PCA, sNMF admixture (K-scan), UMAP, DAPC, cross-method consensus
  (ARI/NMI), bootstrap-supported NJ trees, Patterson f3 admixture tests, anchored to the
  3,000 Rice Genomes reference (XI/GJ/cA/cB subpopulations).
- **Diversity & mating system**: He/Ho/Fis, selfing rate, rarefied allelic richness,
  private alleles, π/θ_W/Tajima's D on a no-MAF-filter site set (the MAF-filter artifact
  is documented).
- **Genome history**: ROH/F_ROH (with a residually heterozygous group isolated), LD-based
  Ne (selfing-adjusted), folded SFS + Stairway Plot demography (Set 1).
- **Differentiation & selection**: Weir–Cockerham θ with block-bootstrap CIs, Nei Gst,
  one- and two-level AMOVA, LD decay, pcadapt + phased iHS/XP-EHH scans with
  convergence-based interpretation.
- **Genomic↔phenomic concordance**: does image-based phenomics (204 canopy features)
  recover the genetic population structure? Mantel/partial Mantel, Procrustes/PROTEST,
  supervised classification vs permutation nulls, MMRR, feature-family attribution,
  size-leakage and residualisation tests.
- **Breeding resource**: the narrow base quantified in Ne units, plus a
  diversity-maximising core collection.

Framing rule used throughout: the two panels are disjoint and on different platforms;
each is analysed separately and the results compared; only the 3K-RGP-common-SNP
co-analysis merges data (with stated caveats). They are never presented as "one
300-accession panel on two platforms."

## Repository layout

```text
README.md              this file
MANIFEST.md            FLOAT → SCRIPT → INPUTS → OUTPUTS provenance map for every figure/table
run_all.sh             master script (one command; stages [0]–[26] + audit gates)
environment.yml        pinned R + Python environment (provenance note inside)
restore_env.sh         rebuild the pinned environment (conda)
CITATION.cff, LICENSE
analysis/
  scripts/             numbered pipeline 01→38 + gates 90–92 (Python + R); paths.py/_paths.R
  results/
    figures/{main,supp}/   Fig01–Fig10 + SuppFig01–05 (PNG + PDF)
    tables/                result tables produced by the stages
    supp_tables/           SuppTable_S1–S20 (CSV)
    source_data/           Source_Data.xlsx (one sheet per data-bearing main figure)
```

## Reproducibility

**What's included vs. what isn't.** Every derived result table and figure file produced by this
pipeline is committed here (and archived as the versioned data asset of release v1.0.0), so the
full result set can be inspected without running anything. The raw genotype/phenotype inputs and
the PLINK 1.9 binary that `run_all.sh` needs to *regenerate* those results from scratch are
restricted/institutional sources, not distributed in this repository — see the paper's Data
Availability statement, and `analysis/data/external/README.md` for exactly which files/env vars
each stage expects. Running `run_all.sh` without them will stop with a clear
`[external data missing]` error naming the missing input, rather than a bare crash.

```bash
bash restore_env.sh      # once: rebuild the pinned environment
bash run_all.sh          # regenerate every figure/table (needs the external data above)
```

- `MANIFEST.md` maps every float to its producing script and inputs, and lists the external-data
  dependencies explicitly.
- The inventory audit (`analysis/scripts/12_audit.py`) checks the full expected figure/table
  inventory and exits non-zero on failure.
- External tools (PLINK 1.9, Beagle 5.5, Stairway Plot 2, Java 17) are resolved through
  `analysis/scripts/paths.py` / `_paths.R`; R package requirements are listed in
  `environment.yml`.

## Citation and license

Please cite the associated paper (details on publication) and this archive; see
`CITATION.cff`. Code is released under the MIT License.
