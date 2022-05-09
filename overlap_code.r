library(GenomicRanges)
library(rtracklayer)
library(methylKit)
library(genomation)
library(Biostrings)
# library(stringr)

# Input: 2 GRange Objects
# Output: 1 GRange Object (same as first input with additional Column)
# Parameter
#  naming: Name of the new Column
overlap_annot <- function(interval, annot, naming, ColName="Name") {
  interval@elementMetadata[[naming]] <- rep("none", length(interval))
  overlap <- findOverlaps(annot,interval)
  uniqueannot <- !duplicated(queryHits(overlap))
  interval@elementMetadata[[naming]][subjectHits(overlap)[uniqueannot]] <- as.character(annot[queryHits(overlap)[uniqueannot],]@elementMetadata[[ColName]])
  return(interval)
}

PromDiff <- calculateDiffMeth(DMpromoters, overdispersion="MN", mc.cores = 6)

PromDiff_hyper <- as(getMethylDiff(PromDiff,difference=mydiffc,qvalue=mydiffp,type="hyper"), "GRanges")
PromDiff_hyper <- overlap_annot(PromDiff_hyper, Promoters,"ID", "ID")
