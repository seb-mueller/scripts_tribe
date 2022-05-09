source(file.path(path_base, "scripts/functions.R"))
source(file.path(path_base, "scripts/R-functions-seb.R"))
library(GenomicRanges)
library(datasets)
library(rtracklayer)
library(methylKit)
# library(genomation)
library(stringr)
library(Biostrings)
# library(topGO)
library(readr)
library(ggplot2)
library(purrr)
library(forcats)
library(viridis)
library(magrittr)
library(ggpubr)
library(purrr)
library(rlang)
library(tibble)
library(tidyr)
library(RColorBrewer)


gglayers2 <- list(
                  theme_bw(),
                  scale_y_continuous(breaks = scales::pretty_breaks(n = 8)),
                  theme(
                        panel.grid.minor.y = element_blank(),
                        legend.position = "right",
                        axis.ticks = element_blank(),
                        axis.text.x = element_text(angle = 300, hjust = 0)))

gglayers3 <- append(gglayers2, scale_y_continuous(labels = scales::percent_format(acc = 1, scale = 1)))

mygenomes <- c("SolLyc", "SolPen")
mygenome <- mygenomes[1]
# load(file = file.path(path_comparison, paste0(mygenome, "_bins_plus_statsv16.rdata"))) # slydf
F4s_all <- c("P1512", "P1515","P2561","P2562","P3611","P3612","P4041","P4042")
F4s_rnaseq <- F4s_all[! F4s_all %in% "P1515"]
F4s_meth <-  F4s_all[! F4s_all %in% "P2562"]
# since we don't have genotype info for P1515:
F4s_meth2 <-  F4s_all[! F4s_all %in% c("P2562", "P1515")]

# F4_shapes <- 1:8
# names(F4_shapes) <- F4s_all
# F4_shapes <- structure(1:8, .Names = c("P1512", "P1515", "P2561", "P2562", "P3611", "P3612", "P4041", "P4042"))
# TRIBE_colors <- brewer.pal(5,"Set3")
# names(TRIBE_colors) <- c("non-DESL", "inbetween","none","+DESL", "-DESL")
# TRIBE_colors["inbetween"] <- "#bebada"
# TRIBE_colors["none"] <- "#8dd3c7"
# TRIBE_colors["non-DESL"] <-"#ffffb3"
# TRIBE_colors["non-DEG"] <- "#ffffb3"
# TRIBE_colors["non-DET"] <- "#ffffb3"
# TRIBE_colors["non-DMR"] <- "#ffffb3"
# TRIBE_colors["+DEG"] <- TRIBE_colors["+DESL"]
# TRIBE_colors["-DEG"] <- TRIBE_colors["-DESL"]
# TRIBE_colors["+DET"] <- TRIBE_colors["+DESL"]
# TRIBE_colors["-DET"] <- TRIBE_colors["-DESL"]
# TRIBE_colors["+DMR"] <- TRIBE_colors["+DESL"]
# TRIBE_colors["-DMR"] <- TRIBE_colors["-DESL"]
# dput(TRIBE_colors)
# dput(F4_shapes)

# easy import:
TRIBE_colors <- c(`non-DESL` = "#ffffb3", inbetween = "#bebada", none = "#8dd3c7",
`+DESL` = "#FB8072", `-DESL` = "#80B1D3", `non-DEG` = "#ffffb3",
`non-DET` = "#ffffb3", `non-DMR` = "#ffffb3", `+DEG` = "#FB8072",
`-DEG` = "#80B1D3", `+DET` = "#FB8072", `-DET` = "#80B1D3", `+DMR` = "#FB8072",
`-DMR` = "#80B1D3")
F4_shapes <- structure(1:8, .Names = c("P1512", "P1515", "P2561", "P2562",
"P3611", "P3612", "P4041", "P4042"))

mybreaks <- c(0, 0.05, 0.9, 1)
mygenomes <- c("SolLyc", "SolPen")
intervaltype <- "bins"
# only Sly for now:
mygenome <- mygenomes[1]
# attach(.sebenv) 
path_scripts <- file.path(path_base, "scripts")
path_comparison <- file.path(path_base, "comparisons_using_big_table/")
path_meth <- file.path(path_base, "bsseq_2018/analysis_merged_Sly_Spen")
path_srnas_bins <- file.path(path_base, "srnas/analysis_bins/")
path_srnas<- file.path(path_base, "srnas")
path_base_rnaseq <-  "/projects/TRIBE/RNA_seq/"
path_map  <- file.path(path_base_rnaseq, "20170901_RNAseq_parents_and_F4_preprocessed_by_Sara/mapped_data/STAR_SL3.0_plus_organelles")
myF4s <-   c("P1512", "P2562", "P2561","P3611","P3612","P4041", "P4042")
F4s_all <- c("P1512", "P1515","P2561","P2562","P3611","P3612","P4041","P4042")
F4s_srnas <- F4s_all[! F4s_all %in% "P1515"]
F4s_meth <-  F4s_all[! F4s_all %in% "P2562"]
F4s_rnaseq <- F4s_all[! F4s_all %in% "P1515"]
contexts  <-  c("CpG", "CHH", "CHG")
gtfilter <- "00"

genes_sly <- "/data/public_data/tomato/ITAG3.2/ITAG3.2_gene_models.gtf2"
genes_sly_fasta <- "/data/public_data/tomato/ITAG3.2/ITAG3.2_CDS.fasta"
TEs_path <- "/data/public_data/tomato/ITAG3.2/ITAG3.2_REPET_repeats_agressive.gtf"
genome_sly_path <- "/data/public_data/tomato/S_lycopersicum/S_lycopersicum_chromosomes.3.00.fa"
genome_penn_path <- "/data/public_data/tomato/Solanum_pennellii_complete/Solanum_pennellii_genome/Solanum_pennellii_genome.fasta"
# tomato_fa <- readDNAStringSet(genome_sly_path, format="fasta")
# names(tomato_fa) <-  str_replace(names(tomato_fa)," ","")

# TODO into the annotation
mirRNAs_path <- "/data/public_data/tomato/additional_resources/miRNAs_lyc_and_penn/mapping_sebastian/miRNAs_sly_penn_multi.bed"
