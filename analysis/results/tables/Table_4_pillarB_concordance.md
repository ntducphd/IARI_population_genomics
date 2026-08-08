**Table 4a. Mantel and partial-Mantel tests: genomic, phenomic, and trait distances.**

| Panel   |   n | Comparison                      |   Mantel r |     P | Bootstrap 95% CI   |
|:--------|----:|:--------------------------------|-----------:|------:|:-------------------|
| Set1    | 146 | genomic~phenomic                |     0.1597 | 0.001 | [0.122, 0.203]     |
| Set1    | 146 | genomic~trait                   |     0.0617 | 0.001 | --                 |
| Set1    | 146 | phenomic~trait                  |     0.25   | 0.001 | --                 |
| Set1    | 146 | genomic~phenomic|trait(partial) |     0.1493 | 0.001 | --                 |
| Set1    | 146 | genomic~trait|phenomic(partial) |     0.0228 | 0.135 | --                 |
| Set2    | 147 | genomic~phenomic                |     0.0867 | 0.004 | [0.029, 0.183]     |
| Set2    | 147 | genomic~trait                   |     0.0915 | 0.001 | --                 |
| Set2    | 147 | phenomic~trait                  |     0.1353 | 0.014 | --                 |
| Set2    | 147 | genomic~phenomic|trait(partial) |     0.0753 | 0.012 | --                 |
| Set2    | 147 | genomic~trait|phenomic(partial) |     0.0809 | 0.002 | --                 |

**Table 4b. Procrustes/PROTEST: genomic-PCA vs phenomic-PCA.**

| Panel   |   n |   Procrustes M2 |   Permutation P |   N permutations |
|:--------|----:|----------------:|----------------:|-----------------:|
| Set1    | 146 |          0.8047 |           0.001 |              999 |
| Set2    | 147 |          0.8564 |           0.001 |              999 |

**Table 4c. Supervised classification: phenomic features -> admixture cluster.**

| Panel   |   n |   Admixture clusters (K) |   CV folds |   RF CV accuracy |   Majority-class baseline |
|:--------|----:|-------------------------:|-----------:|-----------------:|--------------------------:|
| Set1    | 146 |                        7 |          2 |           0.4521 |                    0.226  |
| Set2    | 147 |                        9 |          5 |           0.2993 |                    0.1905 |

**Table 4d. Feature-type attribution: which imaging channel carries genomic signal.**

| Panel   | Feature family   |   N features |   Genomic Mantel r |     P |
|:--------|:-----------------|-------------:|-------------------:|------:|
| Set1    | Size_morphology  |          111 |             0.1587 | 0.001 |
| Set1    | Colour           |           81 |             0.1262 | 0.001 |
| Set1    | NIR              |           12 |             0.091  | 0.001 |
| Set2    | Size_morphology  |          111 |             0.1396 | 0.001 |
| Set2    | Colour           |           81 |             0.0269 | 0.348 |
| Set2    | NIR              |           12 |             0.1442 | 0.001 |

**Table 4e. Structure-as-confounder test: phenomic~trait, raw vs partialled for genomic distance.**

| Panel   |   Raw phenomic~trait r |   Raw P |   Partial r (| genomic) |   Partial P |   % change |
|:--------|-----------------------:|--------:|------------------------:|------------:|-----------:|
| Set1    |                 0.25   |   0.001 |                  0.2437 |       0.001 |        2.5 |
| Set2    |                 0.1353 |   0.014 |                  0.1284 |       0.02  |        5.1 |
