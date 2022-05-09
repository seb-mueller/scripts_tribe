#' ---
#' title:  Common objects and snippets for tomato analysis
#' author: Sebastian Mueller (sm934)
#' date:   2019-04-25
#' ---
# .libPaths(.libPaths()[3:2])
library(GenomicAlignments)
library(grid)
library(rtracklayer)
library(RColorBrewer)
library(edgeR)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(forcats)
library(readr)
library(readxl)

overlap_annot <- function(query, subject, naming, ColName="Name", uniqueSubject = TRUE) {
  # query=head(tmp, 2)
  # subject=TEs[["Slyc"]]
  # naming="test"
  # ColName="Name"
  query@elementMetadata[[naming]] <- rep("none", length(query))
  overlap <- findOverlaps(subject,query)
  # overlap <- head(overlap)
  if (uniqueSubject) {
    uniquesubject <- !duplicated(queryHits(overlap))
    query@elementMetadata[[naming]][subjectHits(overlap)[uniquesubject]] <-
      as.character(subject[queryHits(overlap)[uniquesubject],]@elementMetadata[[ColName]])
  } else {
    query@elementMetadata[[naming]][subjectHits(overlap)] <-
      as.character(subject[queryHits(overlap),]@elementMetadata[[ColName]])
  }
  return(query)
}
# > tmp[1,]
# GRanges object with 1 range and 4 metadata columns:
#        seqnames    ranges strand |          ID       genome        type
#           <Rle> <IRanges>  <Rle> | <character>     <factor> <character>
#   [1] SL3.0ch00     1-200      * | SL3.0ch00-1 SL3.0-Genome         200
#                                                   TEs_name
#                                                <character>
#   [1] ms301134_SL3.0ch00_RLX-incomp_SL3_120k-B-R4111-Map11
#   -------
#   seqinfo: 28 sequences from an unspecified genome
# > sort(TEs[["Slyc"]])[1,]
# GRanges object with 1 range and 12 metadata columns:
#        seqnames    ranges strand |   source     type     score     phase
#           <Rle> <IRanges>  <Rle> | <factor> <factor> <numeric> <integer>
#   [1] SL3.0ch00 1989-5893      + |    REPET    match      <NA>      <NA>
#       Superfamily                                             Target
#       <character>                                        <character>
#   [1]       Gypsy RLX-incomp_SL3_120k-B-R232-Map6_reversed 1175 5056
#                                                                ID       Order
#                                                       <character> <character>

# computes average methylation for each sample in a "methylBase" object:
# by Sebastian Mueller
# e.g. meth:
# methylBase object with 6 rows
# --------------
#   chr start end strand coverage1 numCs1 numTs1 coverage2 numCs2 numTs2 coverage3 numCs3 ...
# 1   1   110 110      *         3      3      0        26     25      1         2      2
# 2   1   115 115      *         4      4      0        50     49      1         3      3
#
# usage: meth_summary(meth)
# output: 
# sample1   sample2 
# 0.7661584 0.7459081
meth_summary <- function(mymeth) {
  require(dplyr)
  mynames <- mymeth@sample.ids
  mysummary <- as.data.frame(mymeth) %>%
    select(starts_with("num")) %>%
    colSums(na.rm = TRUE) %>%
    tapply(rep(seq_along(mynames), each = 2), function(x) x[1] / sum(x))
  names(mysummary) <- mynames
  return(mysummary)
}

# setwd("/data/public_data/tomato/additional_resources/R_code_and_objects/")
# load(file = "R-objects-grafting.rdata")

myPaths <- vector()
myPaths["genome"] <- "/data/public_data/tomato/merged_genomes/genomes_SL30_penn_organelles_merged.fa"

myPaths["genome_penn"] <- "/data/public_data/tomato/Solanum_pennellii_complete/Solanum_pennellii_genome/Solanum_pennellii_genome.fasta"

myPaths["genome_Slyc"] <- "/data/public_data/tomato/S_lycopersicum/S_lycopersicum_chromosomes.3.00.fa"
# myPaths["genes_Slyc"] <- "/data/public_data/tomato/S_lycopersicum/transcriptome_fasta/Slyc_ITAG3.0.gff"
myPaths["genes_Slyc"] <- "/data/public_data/tomato/ITAG3.2/ITAG3.2_gene_models.gff_sorted.gff"
myPaths["genes_Slyc_gtf"] <- "/data/public_data/tomato/ITAG3.2/ITAG3.2_gene_models.gtf2"
myPaths["genes_Slyc_miRNAs"] <- "/data/public_data/tomato/mirbase/lift_over_to_tomato_assembly_3.0"
# myPaths["genes_Spen"] <- "/data/public_data/tomato/Solanum_pennellii_complete/transcriptome_fasta/Spen_gene_models_10.gff"
myPaths["genes_Spen"] <- "/data/public_data/tomato/Solanum_pennellii_complete/Solanum_pennellii_annotation/spenn_v20_gene_models_annot_sorted.gff"
myPaths["genes_Spen_gtf"] <- "/data/public_data/tomato/Solanum_pennellii_complete/Solanum_pennellii_annotation/spenn_v20_gene_models_annot_sorted.gtf"
myPaths["TEsSlyc"] <- "/data/public_data/tomato/ITAG3.2/ITAG3.2_REPET_repeats_agressive.gff"
myPaths["miRNAs_Slyc"]  <- "/home/seb/workspace/tribe/tomato/sl30b_sly_mod.gff3"
# myPaths["TEsSlyc"] <- "/data/public_data/tomato/ITAG3.2/ITAG3.2_RepeatModeler_repeats_light.gff"

Genes_Slyc <- import.gff3(myPaths["genes_Slyc"])
Genes_Spen <- import.gff3(myPaths["genes_Spen"])
miRNA_Slyc <- import.gff3(myPaths["miRNAs_Slyc"])
table(Genes_Spen$type)
table(Genes_Slyc$type)

# # filter out useless colums
# elementMetadata(Genes_Spen) <- 
#   elementMetadata(Genes_Spen) %>%
#   as.data.frame() %>%
#   select(-contains("_Note"))

GOlist <- TEs <- Promoters <- genes <- miRNAs <-  list()
TEs_Slyc <- import.gff3(myPaths["TEsSlyc"])
TEs[["Slyc"]] <- sort(TEs_Slyc)
miRNAs[["Slyc"]] <- miRNA_Slyc
## no TEs for Penn to be found! Running repeatmasker ourselves (or REPET)?
genes[["Slyc"]][32831]$Note <- " " # only one gene empty: fixing
genes[["Slyc"]] <- Genes_Slyc[Genes_Slyc@elementMetadata$type == "mRNA", ]

genes[["Slyc"]]$Note2 <- sapply(genes[["Slyc"]]$Note, function(x) x[[1]])
genes[["Slyc"]]$Note2 <- as.character(genes[["Slyc"]]$Note)
genes[["Slyc"]]$Ontology_term2 <- sapply(genes[["Slyc"]]$Ontology_term, function(x) paste0(x, collapse = ";"))
genes[["Slyc"]]$Ontology_term3 <- sapply(genes[["Slyc"]]$Ontology_term, function(x) x[1])

GOlist[["Slyc"]] <- genes[["Slyc"]]$Ontology_term3
names(GOlist[["Slyc"]]) <- genes[["Slyc"]]$Name

genes[["Spen"]] <- Genes_Spen[Genes_Spen@elementMetadata$type=="gene", ]
Promoters[["Slyc"]] <- promoters(genes[[ "Slyc" ]])
Promoters[["Spen"]] <- promoters(genes[[ "Spen" ]])

tomato_merged <- readDNAStringSet(myPaths["genome"])
chr_sizes <- width(tomato_merged)
names(chr_sizes) <- names(tomato_merged)
chr_sizes
#     SL3.0ch00    SL3.0ch01     SL3.0ch02     SL3.0ch03     SL3.0ch04
#      20852292      98455869      55977580      72290146      66557038
# ...
#    Spenn-ch12 mitochondrion   chloroplast
#      83305730        586297        155461

#------------------------------------ adding heter/eu-chromatin information from cmt3 paper
euchr <- read.csv(file = "/data//public_data/tomato/additional_resources/heterochromatin_definition_cmt3_paper/SL3_euchromatin_chipseq_ratio.csv")
hetchr <- read.csv(file = "/data//public_data/tomato/additional_resources/heterochromatin_definition_cmt3_paper/SL3_heterochromatin_chipseq_ratio.csv")
euchr_gr <- makeGRangesFromDataFrame(euchr, keep.extra.columns=TRUE)
hetchr_gr <- makeGRangesFromDataFrame(hetchr, keep.extra.columns=TRUE)
euchr_gr$state <- "Euchromatin"
hetchr_gr$state <- "Heterochromatin"
chromatin <- c(euchr_gr, hetchr_gr)

#------------------------------------ define bins
# # binning genome into various sizes
bins <- list("2k" = list(width = 2e3), 
             "200" = list(width = 2e2), 
             "500" = list(width = 5e2), 
             "50k" = list(width = 5e5))
# mywidth <- "200"
  tmp <- bins[[mywidth]][["gr"]]
for (mywidth in names(bins)) {
  # mywidth <- names(bins)[2]
  tmp <- tileGenome(seqlengths =  chr_sizes,
                    tilewidth  =  bins[[mywidth]][["width"]],
                    cut.last.tile.in.chrom = TRUE)
  tmp$ID <- paste(seqnames(tmp),start(tmp), sep= "-")
  tmp$genome <- as.factor(paste(substr(seqnames(tmp),1,5),"Genome",sep="-"))
  # tmp$genome <- reshape::combine_factor(tmp$genome,c(1,1,2,3))
  # levels(tmp$genome)[1] <- "both"
  tmp$type <- mywidth
  # tmp <- overlap_annot(tmp, TEs[["Slyc"]], "TEs_name", "Name")
  tmp <- overlap_annot(tmp, chromatin, "chromatin_state", "state", unique = FALSE)
  tmp <- overlap_annot(tmp, TEs[["Slyc"]], "TEs_name", "Name", unique = FALSE)
  tmp <- overlap_annot(tmp, TEs[["Slyc"]], "TEs_Order", "Order", unique = FALSE)
  tmp <- overlap_annot(tmp, TEs[["Slyc"]], "TEs_Class", "Class", unique = FALSE)
  tmp <- overlap_annot(tmp, TEs[["Slyc"]], "TEs_Superfamily", "Superfamily", unique = FALSE)
  tmp <- overlap_annot(tmp, genes[["Slyc"]], "Genes_Slyc", "Name", unique = FALSE)
  tmp <- overlap_annot(tmp, genes[["Spen"]], "Genes_Spen", "Name", unique = FALSE)
  tmp <- overlap_annot(tmp, genes[["Slyc"]], "Ontology_term_Slyc", "Ontology_term2", unique = FALSE)
  tmp <- overlap_annot(tmp, genes[["Slyc"]], "Note", "Note2", unique = FALSE)
  tmp <- overlap_annot(tmp, Promoters[["Slyc"]], "Promoters_Slyc", "Name", unique = FALSE)
  tmp <- overlap_annot(tmp, Promoters[["Spen"]], "Promoters_Spen", "Name", unique = FALSE)
  tmp <- overlap_annot(tmp, miRNAs[["Slyc"]], "miRNAs_Slyc", "Name", unique = FALSE)
  bins[[mywidth]][["gr"]] <- tmp
}

# save(GOlist, genes, Promoters, TEs, tomato_merged, bins, chr_sizes, myPaths, file = "/data/public_data/tomato/additional_resources/R_code_and_objects/R-objects-bins-SL30-Penn_v6_annot.rdata")
# path_comparison <- "/home/seb/workspace/tribe/comparisons_using_big_table/"
save(GOlist, miRNAs, genes, Promoters, TEs, tomato_merged, bins, chr_sizes, myPaths, file = file.path(path_comparison, paste0("R-objects-bins-SL30-Penn_v7_annot.rdata")))
