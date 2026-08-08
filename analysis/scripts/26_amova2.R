#!/usr/bin/env Rscript
# 26_amova2.R -- [stage 26] two-level hierarchical AMOVA for Set 1:
# among macro-groups (Xian/indica vs Geng/japonica vs circum-Aus vs admixed) /
# among 3K-RGP subpopulations within macro-group / within subpopulation.
#
# Motivation: the stage-05b AMOVA has one level
# (admixture cluster). A second, EXTERNAL hierarchy exists for Set 1 -- the authoritative
# 3,000 Rice Genomes/SNP-Seek subpopulation labels (stage 08a) nest naturally into the classical
# indica/japonica/aus macro-groups -- so a two-level variance partition is computable without
# guessing any provenance. (A landrace/cultivar/BAM provenance AMOVA was considered and dropped:
# no per-accession provenance table exists in the data lake, and inferring classes from ID
# prefixes would violate the no-hand-guessing rule. Likewise isolation-by-distance was dropped:
# 147/150 Set 1 accessions share one country of origin (India), so there is no geographic
# variance to test. Both exclusions are stated in the manuscript.)
#
# Method: pegas::amova on the 1-IBS distance, nested formula dist ~ macro/subpop, 999 perms.
# Groups with n < 3 are pooled into "other" (variance components unstable otherwise; reported).
#
# Outputs (analysis/results/tables/):
#   amova2_set1.csv — variance component, df, sigma2, percent, p (where testable)
# Run via PowerShell (house rule).

suppressMessages({ library(pegas) })

args      <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPTS   <- dirname(normalizePath(this_file))
ROOT      <- normalizePath(file.path(SCRIPTS, "..", ".."))
INTERIM   <- file.path(ROOT, "analysis", "data", "interim")
TAB       <- file.path(ROOT, "analysis", "results", "tables")

ibs <- readRDS(file.path(INTERIM, "ibs_dist_set1.rds"))
d   <- as.dist(ibs)
ids <- attr(ibs, "Labels"); if (is.null(ids)) ids <- rownames(as.matrix(ibs))

sub <- read.csv(file.path(TAB, "subpop_assignment_set1.csv"))
lab <- setNames(as.character(sub$subpopulation), as.character(sub$sample))
subpop <- lab[ids]

macro_of <- function(s) {
  if (is.na(s)) return(NA_character_)
  if (grepl("^(indx|ind[0-9]|XI)", s, ignore.case = TRUE)) return("Indica")
  if (grepl("^(temp|trop|subtrop|japx|GJ)", s, ignore.case = TRUE)) return("Japonica")
  if (grepl("^(aus|cA)", s, ignore.case = TRUE)) return("Aus")
  if (grepl("admix", s, ignore.case = TRUE)) return("Admixed")
  "other"
}
macro <- vapply(subpop, macro_of, character(1))

keep <- !is.na(subpop) & !is.na(macro)
# pool subpops with n<3
tab_sub <- table(subpop[keep])
subpop[keep][subpop[keep] %in% names(tab_sub[tab_sub < 3])] <- "pooled_small"

dm <- as.matrix(d)[keep, keep]
d2 <- as.dist(dm)
g1 <- factor(macro[keep]); g2 <- factor(paste(macro[keep], subpop[keep], sep = ":"))

cat(sprintf("[26] Set1 n=%d; macro groups: %s\n", sum(keep),
            paste(sprintf("%s(%d)", levels(g1), table(g1)), collapse = ", ")))

res <- pegas::amova(d2 ~ g1/g2, nperm = 999)
tab <- res$tab
sig2 <- setNames(res$varcomp$sigma2, rownames(res$varcomp))
pct  <- 100 * sig2 / sum(sig2)
pv   <- res$varcomp$P.value

out <- data.frame(component = rownames(res$varcomp), sigma2 = sig2,
                  percent = pct, p = pv)
write.csv(out, file.path(TAB, "amova2_set1.csv"), row.names = FALSE)
print(out)
cat("[26] done -> amova2_set1.csv\n")
