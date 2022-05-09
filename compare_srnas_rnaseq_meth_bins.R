#------------------------------------ done once:
# subsetting bins to only look at Slyc genome
# from /data/public_data/tomato/additional_resources/R_code_and_objects/bins_tomato.r
load("/data/public_data/tomato/additional_resources/R_code_and_objects/R-objects-bins-SL30-Penn_v6_annot.rdata")
# load("/home/sm934/hydro_data/public_data/tomato/additional_resources/R_code_and_objects/R-objects-bins-SL30-Penn_v3_annot.rdata")

sly_bins <- subset(bins[["200"]]$gr, genome == "SL3.0-Genome")
# save(sly_bins, file = file.path(path_comparison, paste0(mygenome, "_bins_lyc.rdata")))

#------------------------------------ start here 
path_base <-  "/projects/TRIBE"
# laptop
path_base <-  "/home/seb/workspace/tribe"
path_base_rnaseq <-  "/projects/TRIBE/RNA_seq/"
source(file.path(path_base, "scripts/project_settings.R"))

setwd(file.path(path_comparison))

# load(file = paste0(mygenome, "_bins_lyc.rdata"))
load(file = file.path(path_comparison, paste0(mygenome, "_bins_lyc.rdata")), 
     envir = binenv <- new.env())

ls(binenv) # [1] "sly_bins"
attach(binenv)

# Note P1511 died and was replaced by P1515

#------------------------------------ adding genotype to bins
hetmaps_condensed2 <- readRDS(file = file.path(path_base, "F4_genotyping", "hetmaps_condensed2.rds"))

libs <- names(hetmaps_condensed2[[mygenome]])
for (i in libs) {
  sly_bins <- overlap_annot(sly_bins, hetmaps_condensed2[[mygenome]][[i]], 
                            paste0("genotype_", i), "score", unique = FALSE)
  # genes[["Spen"]]
}
# sly_bins[2853:2859,]

# adding Hajks genome as annotation
LTR_path <- "/data/public_data/tomato/additional_resources/TE_annotation_Hajk/Slycopersicum_v3_annotation_all.gff"
LTR_path_active <- "/data/public_data/tomato/additional_resources/TE_annotation_Hajk/Slycopersicum_v3_annotation_potentiallyActive_98.gff"
LTR_path_active <- "~/Slycopersicum_v3_annotation_all.gff"

LTRs_sly <- import.gff3(LTR_path_active)

TEs <- TEs[["Slyc"]]
ol <- findOverlaps(TEs, LTRs_sly, minoverlap=50)
tmp1 <- table(TEs[queryHits(ol),]$Superfamily)
tmpsub <- LTRs_sly[subjectHits(ol),]
tmp2 <- table(TEs$Superfamily)

tmp <- as.data.frame(tmp1)
colnames(tmp) <- c("Var1", "Freq2")

left_join(tmp,as.data.frame(tmp2)) %>%
  mutate(ratio = Freq/Freq2)

                                                 # Var1 Freq2   Freq     ratio
# 1                                               CACTA    16   4598 287.37500
# 2                                         Confused_TE   207  43870 211.93237
# 3             Confused_TE_withPotentialHostGeneDomain    63  17032 270.34921
# 4                                               Copia  2368  66977  28.28421
# 5                   Copia_withPotentialHostGeneDomain    64   6979 109.04688
# 6                                                EPRV    78   5699  73.06410
# 7                                               Gypsy  2363 226853  96.00212
# 8                   Gypsy_withPotentialHostGeneDomain    68  14857 218.48529
# 9                                           Harbinger    10   1632 163.20000
# 10                                                hAT    30   5878 195.93333
# 11                    hAT_withPotentialHostGeneDomain     2    771 385.50000
# 12                                           Helitron     7   1703 243.28571
# 13               Helitron_withPotentialHostGeneDomain     1    811 811.00000
# 14                                           HostGene    28  11528 411.71429
# 15                                               LINE    57  18113 317.77193
# 16                   LINE_withPotentialHostGeneDomain     6   1218 203.00000
# 17                                            Mariner     3    953 317.66667
# 18                                               MuDR    31  10343 333.64516
# 19                             PutativeNonAutoClassII    46  17155 372.93478
# 20 PutativeNonAutoClassII_withPotentialHostGeneDomain     7   1429 204.14286
# 21                                        putNA_CACTA     1    288 288.00000
# 22                                          putNA_hAT    14   3557 254.07143
# 23                                    Retrotransposon     2   1108 554.00000
# 24                                                SAT    18   5309 294.94444
# 25                                                SSR    34   6413 188.61765
# 26                                       TIR_MITEov10    49  13938 284.44898
# 27                                          TRIM_LARD     7   1031 147.28571
# 28                                       Unclassified   472 109206 231.36864

# see read_psl.R how the generate those below


#------------------------------------  adding EPRV domain BLAT annotation
path_giri <- file.path(path_base, "EPRV_annotation/EPRV_giri")
load(file = file.path(path_giri, "eprv_domains_blat.rdata")) #eprv_domains_blat

# eprv_domains_blat <- mypsl[width(mypsl) < 1e5,]
TEstmp <- overlap_annot(sly_bins, eprv_domains_blat[eprv_domains_blat$qName == "env",], "domain_env", "matches")
TEstmp <- overlap_annot(TEstmp, eprv_domains_blat[eprv_domains_blat$qName == "inclusion_body",], "domain_inclusion_body", "matches")
TEstmp <- overlap_annot(TEstmp, eprv_domains_blat[eprv_domains_blat$qName == "movement_protein",], "domain_movement_protein", "matches")
sly_bins <- overlap_annot(TEstmp, eprv_domains_blat[eprv_domains_blat$qName == "pol",], "domain_pol", "matches", unique = FALSE) # 406

with(as.data.frame(sly_bins), table(domain_pol != "none", domain_env != "none"))
with(as.data.frame(sly_bins), table(domain_pol != "none", domain_inclusion_body != "none"))
with(as.data.frame(sly_bins), table(domain_pol != "none", domain_movement_protein != "none"))
with(as.data.frame(sly_bins), table(domain_env != "none", domain_movement_protein != "none"))

with(as.data.frame(sly_bins), table(domain_env != "none", srnas_unique_P4042_class == "DESL")) #37
with(as.data.frame(sly_bins), table(domain_pol != "none", grnas_unique_P4042_class == "DESL")) #15
with(as.data.frame(sly_bins), table(domain_inclusion_body != "none", srnas_unique_P4042_class == "DESL")) #47
with(slydf, table(domain_movement_protein != "none", srnas_unique_P4042_class == "DESL")) #14


#------------------------------------ srnas stats to bins
# srna_cnts_bin <- list(mygenomes[1], mygenomes[2])
srna_prefix <- "sRNA_bins200_"
myparent <- "PenneC"
# myparent <- "M82C"
load(file = file.path(path_srna, paste0(srna_prefix, mygenome,myparent, "_bins_raw_processed.rdata"))) #srna_cnts_bin
load(file = file.path(path_srna, paste0(srna_prefix, mygenome, "_bins_raw_processed_unique.rdata"))) #srna_cnts_unique_bin
load(file = file.path(path_srna, paste0(srna_prefix, mygenome, "_bins_raw_processed_dcl2.rdata")) #lrt_dcl2, cnts_dcl2, 
load(file = file.path(path_srnas, paste0(srna_prefix, mygenome, "_bins_raw_processed_tmv.rdata"))) #lrt_tmv, mycounts_sub 

# load(file = paste0(prefix, mygenome, "_counts_misc.rdata"))
# mycounts, mycounts21, mycounts24, 

for (myF4 in F4s_srnas) {
  # myF4 <-  F4s_srnas[7]
  lrt <- srna_cnts_bin[[myF4]]

  tmp2 <- lrt$table[, c("logFC", "logCPM", "PValue")] %>%
    mutate(FDR = round(p.adjust(PValue, method = "BH"), 5)) %>%
    mutate(class = cut(FDR, breaks = mybreaks,
                       labels = c("DESL", "inbetween", "non-DESL"),
                       include.lowest = TRUE)) %>%
    mutate(updown = ifelse(logFC < 0, "-", "+")) %>%
    mutate(class2 = ifelse(class == "DESL", paste0(updown, as.character(class)),
                           as.character(class))) %>%
    mutate(logCPM = round(logCPM,2)) %>%
    mutate(logFC = round(logFC,2)) %>%
    dplyr::select(-c(updown, PValue, class))

  colnames(tmp2) <- paste("srnas", myF4,myparent, colnames(tmp2), sep = "_")
  # tmp <- sly_bins
  slydf <- cbind(slydf,  tmp2)
}
 with(slydf, table(srnas_P4042_PenneC_class2, srnas_P4042_class2))

## manually inspecting stats!!
 slydf %>%
  filter(genotype_P4042 == "00") %>%
  filter(srnas_P4042_PenneC_class2 == "-DESL") %>%
  filter(srnas_P4042_class2 == "-DESL") %>%
  # dplyr::count(TEs_Order_renamed)
  dplyr::filter(TEs_Order_renamed == "Pararetrovirus") %>%
  # dplyr::count(srnas_tmv_class2)
  dplyr::count(srnas_dcl2_class2)
  # dplyr::count(srnas_P3611_class2)
  # count(TEs_Order_renamed, genotype_P4042)
  # print(n = Inf)

  table(srnas_P4042_PenneC_class2, srnas_P4042_class2))
# with(slydf2, table(srnas_P1512_PenneC_class2))
# with(slydf2, table(srnas_P1512_class2))

# unique
for (myF4 in F4s_srnas) {
  # myF4 <-  F4s_srnas[7]
  lrt <- srna_cnts_unique_bin[[myF4]]

  tmp2 <- lrt$table[, c("logFC", "logCPM", "PValue")] %>%
    mutate(FDR = round(p.adjust(PValue, method = "BH"), 5)) %>%
    mutate(class = cut(FDR, breaks = mybreaks,
                       labels = c("DESL", "inbetween", "non-DESL"),
                       include.lowest = TRUE)) %>%
    mutate(updown = ifelse(logFC < 0, "-", "+")) %>%
    mutate(class2 = ifelse(class == "DESL", paste0(updown, as.character(class)),
                           as.character(class))) %>%
    mutate(logCPM = round(logCPM,2)) %>%
    mutate(logFC = round(logFC,2)) %>%
    dplyr::select(-updown)

  colnames(tmp2) <- paste("srnas_unique", myF4, colnames(tmp2), sep = "_")
  tmp <- sly_bins
  mcols(sly_bins) <- cbind(mcols(sly_bins),  tmp2)
}
 # tmp <- sly_bins
 # sly_bins <- tmp
 # mcols(sly_bins) <- mcols(sly_bins)[, -grep("srna", colnames(mcols(sly_bins)))]
 tmp <- mcols(sly_bins) %>%
   as.data.frame() %>%
   dplyr::select(matches("srnas_unique*.*class2$"))

 tmp2 <- slydf %>%
   dplyr::select(-matches("srnas_unique*.*class2$"))
 slydf <- cbind(tmp2, tmp)
 # mcols(sly_bins) <- mcols(sly_bins)[, -grep("rnaseq", colnames(mcols(sly_bins)))]

# save(sly_bins, file = file.path(path_comparison, paste0(srna_prefix, mygenome, "_bins_plus_statsv2.rdata")))
# sly_bins$myfilter

 #------------------------------------ adding sRNA counts for each F4 and each srna size {{{
load(file = paste0(prefix, mygenome, "_counts_misc_filter.rdata"))
# mycounts21only, mycounts22only, mycounts23only, mycounts24only, 
myF4 ="P4042"
# gtfilter = "00"
myprefix <- "srnas_"
suffix <- "bp_count"

# $lib.sizes from srnas_analysis_tribe_sm934.R
mean_lib_sizes <- mean(meta_sub$lib.sizes)

for (myF4 in c("M82C", "PenneC", F4s_srnas)) {
  idx <- which(meta_sub$plant == myF4)
  meta_sub_F4 <- meta_sub[idx,]
  lib_sizes_F4 <- meta_sub[idx,]$lib.sizes
 # add normalized counts for each size class
  slydf[,paste0(myprefix, myF4, "_20" , suffix)] <- rowMeans(t(t(mycounts20only[, idx] * mean_lib_sizes )/ lib_sizes_F4))
  slydf[,paste0(myprefix, myF4, "_21" , suffix)] <- rowMeans(t(t(mycounts21only[, idx] * mean_lib_sizes )/ lib_sizes_F4))
  slydf[,paste0(myprefix, myF4, "_22" , suffix)] <- rowMeans(t(t(mycounts22only[, idx] * mean_lib_sizes )/ lib_sizes_F4))
  slydf[,paste0(myprefix, myF4, "_23" , suffix)] <- rowMeans(t(t(mycounts23only[, idx] * mean_lib_sizes )/ lib_sizes_F4))
  slydf[,paste0(myprefix, myF4, "_24" , suffix)] <- rowMeans(t(t(mycounts24only[, idx] * mean_lib_sizes )/ lib_sizes_F4))
  slydf[,paste0(myprefix, myF4, "_25" , suffix)] <- rowMeans(t(t(mycounts25only[, idx] * mean_lib_sizes )/ lib_sizes_F4))
  # slydf[,paste0(myprefix, myF4, "_25" , suffix)] <- rowSums(mycounts25only[, which(meta_sub$plant == myF4)])
}

# dcl2 counts
# load(file = paste0(prefix, mygenome, "_counts_dcl2.rdata"))
grep("^wt-",colnames(mycounts20dcl2))
slydf[,paste0("srnas_dcl2_", "d2Mutant", "_24" , suffix)] <- rowMeans(mycounts24dcl2[, 4:5])
slydf[,paste0("srnas_dcl2_", "d2Mutant", "_22" , suffix)] <- rowMeans(mycounts22dcl2[, 4:5])
slydf[,paste0("srnas_dcl2_", "d2Mutant", "_21" , suffix)] <- rowMeans(mycounts21dcl2[, 4:5])
slydf[,paste0("srnas_dcl2_", "d2Mutant", "_25" , suffix)] <- rowMeans(mycounts25dcl2[, 4:5])
slydf[,paste0("srnas_dcl2_", "d2Mutant", "_23" , suffix)] <- rowMeans(mycounts23dcl2[, 4:5])
slydf[,paste0("srnas_dcl2_", "d2Mutant", "_20" , suffix)] <- rowMeans(mycounts20dcl2[, 4:5])
slydf[,paste0("srnas_dcl2_", "WT", "_20" , suffix)] <- rowMeans(mycounts20dcl2[, 1:2])
slydf[,paste0("srnas_dcl2_", "WT", "_23" , suffix)] <- rowMeans(mycounts23dcl2[, 1:2])
slydf[,paste0("srnas_dcl2_", "WT", "_25" , suffix)] <- rowMeans(mycounts25dcl2[, 1:2])
slydf[,paste0("srnas_dcl2_", "WT", "_21" , suffix)] <- rowMeans(mycounts21dcl2[, 1:2])
slydf[,paste0("srnas_dcl2_", "WT", "_22" , suffix)] <- rowMeans(mycounts22dcl2[, 1:2])
slydf[,paste0("srnas_dcl2_", "WT", "_24" , suffix)] <- rowMeans(mycounts24dcl2[, 1:2])
slydf[,paste0("srnas_dcl2_", "tmv", "_20" , suffix)] <- rowMeans(mycounts20dcl2[, 6:7])
slydf[,paste0("srnas_dcl2_", "tmv", "_23" , suffix)] <- rowMeans(mycounts23dcl2[, 6:7])
slydf[,paste0("srnas_dcl2_", "tmv", "_25" , suffix)] <- rowMeans(mycounts25dcl2[, 6:7])
slydf[,paste0("srnas_dcl2_", "tmv", "_21" , suffix)] <- rowMeans(mycounts21dcl2[, 6:7])
slydf[,paste0("srnas_dcl2_", "tmv", "_22" , suffix)] <- rowMeans(mycounts22dcl2[, 6:7])
slydf[,paste0("srnas_dcl2_", "tmv", "_24" , suffix)] <- rowMeans(mycounts24dcl2[, 6:7])



# }}}

#------------------------------------ adding meth stats to bins
# load(file = file.path(path_comparison, paste0(srna_prefix, mygenome, "_bins_plus_statsv2.rdata")), 
     # envir = binenv <- new.env())
# attach(binenv)
# TODO [1] "P2562" is dodgy gor Penn CHH
mygenome  <- "SolLyc"
for (myplant in F4s_meth) {
  for (context in c("CpG", "CHH", "CHG")) {
    # context="CHH"
    # myplant="P4041"
    prefix2 <- paste(context, mygenome, sep = "_")
    message("reading in: ", prefix2)
    prefix3 <- paste(intervaltype, context, mygenome, sep = "_")
    load(file = file.path(path_meth, "bins_Slyc", 
                          paste0(prefix3, "_", myplant, "_stats_interval.rdata")))
    # loading intervalDiff, intervalDiff_allgr
    # intervalDiff contains all bins with a certain coverage (not only DMRs but all)!
    # The idea is to add all info to the big table and filter for DMRs later
    intervalDiff_all2 <- getMethylDiff(intervalDiff,difference=0,qvalue=1)
    intervalDiff_all2_gr <- as(intervalDiff_all2, "GRanges")[,-1] # excluding pvalue (qvalue should suffice)

    tmp <- as.data.frame(elementMetadata(intervalDiff_all2_gr)) %>%
      mutate(meth.diff = round(meth.diff, 2)) %>%
      mutate(class = cut(qvalue, 
                         breaks = mybreaks,
                         labels = c("DMR", "inbetween", "non-DMR"),
                         include.lowest = TRUE)) %>%
    mutate(updown = ifelse(meth.diff < 0, "-", "+")) %>%
    mutate(class2 = ifelse(class == "DMR", paste0(updown, as.character(class)),
                           as.character(class))) %>%
    mutate(class_stringent = ifelse(class == "DMR" & abs(meth.diff) < 15, "inbetween", class2))
    elementMetadata(intervalDiff_all2_gr) <- tmp
    # merging back meth bins into big table
    # note that meth bins should be a grange subset of big table - grange magic required
    sly_bins <- overlap_annot(sly_bins, intervalDiff_all2_gr, 
                              paste("meth", myplant, context, "qvalue", sep = "_"), "qvalue")
    sly_bins <- overlap_annot(sly_bins, intervalDiff_all2_gr, 
                              paste("meth", myplant, context, "meth.diff", sep = "_"), "meth.diff")
    sly_bins <- overlap_annot(sly_bins, intervalDiff_all2_gr, 
                              paste("meth", myplant, context, "class", sep = "_"), "class")
    sly_bins <- overlap_annot(sly_bins, intervalDiff_all2_gr, 
                              paste("meth", myplant, context, "class2", sep = "_"), "class2")
    sly_bins <- overlap_annot(sly_bins, intervalDiff_all2_gr, 
                              paste("meth", myplant, context, "class_stringent", sep = "_"), "class_stringent")
  }
}

# save(sly_bins, file = file.path(path_comparison, paste0(mygenome, "_bins_plus_statsv3.rdata")))

#------------------------------------ adding RNA-Seq
# the strategy here is to mark every bin with wheater it is overlapping with a DE gene or promotor
# load(file = file.path(path_comparison, paste0(mygenome, "_bins_plus_statsv3.rdata")))
# dir_rnaseq <- file.path(path_base, "RNA_seq/20170901_RNAseq_parents_and_F4_sara_analysis_tom/fastq_pooled/RData"
# dir_rnaseq <- file.path(path_base, "RNA_seq/20170901_RNAseq_parents_and_F4_sara_analysis_tom/fastq_pooled/RData"
dir_rnaseq <- file.path(path_base, "20200217_RNAseq_parents_and_F4_seb_analyis_sara_mapped")
# allels <- "Slyc"; other_allels <- "Spen"; mygenome <- "Slyc"
mygenome <- "Slyc"

prefix <- "RNAseq_genes_"
rna_cnts <- readRDS(file = file.path(dir_rnaseq, paste0(prefix, "SolLyc_edgeR_Robject.rds")))
tes_cnts <- readRDS(file = file.path(dir_rnaseq, paste0("RNAseq_TEs_", "SolLyc_edgeR_Robject.rds")))
myparent <- "PenneC"
path_base_rnaseq <-  "/projects/TRIBE/RNA_seq/"
rna_cntspenn <- readRDS(file = file.path(path_base_rnaseq, paste0(prefix, "_vs_", myparent,  "_SolLyc_edgeR_Robject.rds")))
tes_cntspenn <- readRDS(file = file.path(path_base_rnaseq, paste0("RNAseq_TEs_", "_vs_", myparent,  "_SolLyc_edgeR_Robject.rds")))

for (myF4 in F4s_rnaseq) {
  # myF4 <- F4s_rnaseq[1]
  message(myF4)
  rnaseq_fdr <- rna_cnts[[myF4]]$table_all
  gr <- makeGRangesFromDataFrame(rnaseq_fdr, keep.extra.columns=TRUE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("rnaseq", myF4, "FDR", sep = "_"), "FDR", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("rnaseq", myF4, "logFC", sep = "_"), "logFC", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("rnaseq", myF4, "class2", sep = "_"), "class2", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("rnaseq", myF4, "class", sep = "_"), "class", unique = FALSE)


  # associating logFC also with its promotor in an extra column!!
  # this is to facilitate downstream analysis (e.g. correlating bins sitting in promoter with logFC of its gene)
  gr_prom <- promoters(gr)
  sly_bins <- overlap_annot(sly_bins, gr_prom, 
                            paste("rnaseq", myF4, "logFC_promoter", sep = "_"), "logFC", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr_prom, 
                            paste("rnaseq", myF4, "FDR_promoter", sep = "_"), "FDR", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr_prom, 
                            paste("rnaseq", myF4, "class2_promoter", sep = "_"), "class2", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr_prom, 
                            paste("rnaseq", myF4, "class_promoter", sep = "_"), "class", unique = FALSE)
}
for (myF4 in F4s_rnaseq) {
  # for penn
  message(myF4)
  rnaseq_fdr <- rna_cntspenn[[myF4]]$table_all
  gr <- makeGRangesFromDataFrame(rnaseq_fdr, keep.extra.columns=TRUE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("rnaseq", myF4,myparent, "logFC", sep = "_"), "logFC", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("rnaseq", myF4,myparent, "class2", sep = "_"), "class2", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("rnaseq", myF4,myparent, "FDR", sep = "_"), "FDR", unique = FALSE)
}
for (myF4 in F4s_rnaseq) {
  # for penn
  message(myF4)
  rnaseq_fdr <- tes_cntspenn[[myF4]]$table_all
  gr <- makeGRangesFromDataFrame(rnaseq_fdr, keep.extra.columns=TRUE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("TEs_rnaseq", myF4,myparent, "logFC", sep = "_"), "logFC", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("TEs_rnaseq", myF4,myparent, "class2", sep = "_"), "class2", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("TEs_rnaseq", myF4,myparent, "FDR", sep = "_"), "FDR", unique = FALSE)
}

tmp <- mcols(sly_bins) %>%
  as.data.frame() %>%
  dplyr::select(matches("PenneC")) %>%
  mutate(across(!contains("class"), ~ as.numeric(ifelse(.x == "none", NA, .x))))

slydf <- cbind(slydf,  tmp)

for (myF4 in F4s_rnaseq) {
  # myF4 <- F4s_rnaseq[1]
  message(myF4)
  rnaseq_fdr <- tes_cnts[[myF4]]$table_all
  gr <- makeGRangesFromDataFrame(rnaseq_fdr, keep.extra.columns=TRUE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("TEs_rnaseq", myF4, "FDR", sep = "_"), "FDR", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("TEs_rnaseq", myF4, "logFC", sep = "_"), "logFC", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("TEs_rnaseq", myF4, "class2", sep = "_"), "class2", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("TEs_rnaseq", myF4, "class", sep = "_"), "class", unique = FALSE)
}

# add Hayks TEs
  sly_bins <- overlap_annot(sly_bins, LTRs_sly, 
                            "TEs_LTRpred", "ID", unique = FALSE)
# add D2G : Dcl2 rnaseq DEG
 # lrt_dcl2 <- readRDS(file = file.path(path_base_rnaseq, paste0(prefix, "SolLyc_edgeR_Robject_dcl2_deg.rds"))) #lrt_dcl2
lrt_dcl2 <- readRDS(file = file.path(path_base_rnaseq, paste0("RNAseq_dcl2b_nonvirusSolLyc_edgeR_Robject_dcl2_deg.rds"))) #lrt_dcl2

  rnaseq_fdr <- lrt_dcl2$table_all
  gr <- makeGRangesFromDataFrame(rnaseq_fdr, keep.extra.columns=TRUE)

  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("rnaseq_dcl2", "logFC", sep = "_"), "logFC", unique = FALSE)
  gr_prom <- promoters(gr)
  sly_bins <- overlap_annot(sly_bins, gr_prom, 
                            paste("rnaseq_dcl2", "logFC_promoter", sep = "_"), "logFC", unique = FALSE)
  sly_bins <- overlap_annot(sly_bins, gr, 
                            paste("rnaseq_dcl2", "class2", sep = "_"), "class2", unique = FALSE)
  gr_prom <- promoters(gr)
  sly_bins <- overlap_annot(sly_bins, gr_prom, 
                            paste("rnaseq_dcl2", "class2_promoter", sep = "_"), "class2", unique = FALSE)

  # quick inclusion
  slydf$rnaseq_dcl2_class2_promoter <- sly_bins$rnaseq_dcl2_class2_promoter
  slydf$rnaseq_dcl2_class2 <- sly_bins$rnaseq_dcl2_class2
  slydf$rnaseq_dcl2_logFC_promoter <- sly_bins$rnaseq_dcl2_logFC_promoter
  slydf$rnaseq_dcl2_logFC <- sly_bins$rnaseq_dcl2_logFC

slydf <- slydf  %>%
  mutate(rnaseq_dcl2_logFC = as.numeric(ifelse(rnaseq_dcl2_logFC == "none", NA, rnaseq_dcl2_logFC))) %>%
  mutate(rnaseq_dcl2_logFC_promoter = as.numeric(ifelse(rnaseq_dcl2_logFC_promoter == "none", NA, rnaseq_dcl2_logFC_promoter)))

# table(sly_bins$rnaseq_P4042_FDR)
# > rnaseq_fdr[164,]
#      X1                      X2         X3                          X4 gene_models_Slyc_1.M82C.1 gene_models_Slyc_13.M82C.2 gene_models_Slyc_20.P3612.2
# 164 863 mRNA:Solyc01g005300.3.1 SL3.0ch01+ 216877-217040,218535-220489                         0                          4                         646
#     gene_models_Slyc_8.P3612.1     likes {M82C},{P3612} FDR.{M82C},{P3612} FWER.{M82C},{P3612}    logFC         genename  start    end     chrom strand         FDR
# 164                       1404 0.9930063            2>1        0.003188767           0.4079171 8.417853 Solyc01g005300.1 216877 220489 SL3.0ch01      + 0.003188767

#------------------------------------ done adding RNA-Seq

# some postprocessing and cleaning up the table
# converting character columns into numeric

slydf <- slydf  %>%
  mutate(rnaseq_dcl2_logFC = as.numeric(ifelse(rnaseq_dcl2_logFC == "none", NA, rnaseq_dcl2_logFC))) %>%
  mutate(rnaseq_dcl2_logFC_promoter = as.numeric(ifelse(rnaseq_dcl2_logFC_promoter == "none", NA, rnaseq_dcl2_logFC_promoter)))

tmp
 table2(slydf$rnaseq_dcl2_logFC_promoter<0)
 table2(tmp$rnaseq_dcl2_logFC_promoter<0)

tmpdf <- as.data.frame(elementMetadata(sly_bins)) %>% 
  mutate_at(2:8, as.factor) %>%
  # mutate_at(vars(starts_with("rnaseq")), as.numeric) %>%
  mutate_at(vars(ends_with("diff")), as.numeric) %>%
  mutate_at(vars(ends_with("qvalue")), as.numeric) %>%
  mutate_at(vars(ends_with("PValue")), as.numeric) %>%
  mutate_at(vars(ends_with("FDR")), as.numeric) %>%
  mutate_at(vars(ends_with("FDR_promoter")), as.numeric) %>%
  mutate_at(vars(ends_with("CPM")), as.numeric) %>%
  mutate_at(vars(ends_with("logFC")), as.numeric) %>%
  mutate_at(vars(ends_with("logFC_promoter")), as.numeric) %>%
  mutate_at(vars(starts_with("genotype")), as.factor) %>%
  mutate_at(vars(ends_with("class")), as.factor) %>%
  mutate_at(vars(ends_with("class_stringent")), as.factor) %>%
  mutate_at(vars(ends_with("class2")), as.factor) %>%
  mutate_at(vars(matches("TEs_rnaseq.*class*")),
           ~fct_recode(., "-DET" = "-DEG",
                          "+DET" = "+DEG",
                          "non-DET" = "non-DEG",
                          "DET" = "DEG" ))
  
elementMetadata(sly_bins)  <- tmpdf

save(sly_bins, file = file.path(path_comparison, paste0(mygenome, "_bins_plus_statsv23.rdata")))
# load(file = file.path(path_comparison, paste0("SolLyc", "_bins_plus_statsv23.rdata")))
#

slydf <- as.data.frame(sly_bins) %>% 
  mutate(isPromotor = ifelse(Promoters_Slyc != "none", "Promoter", "non-Promoter")) %>%
  mutate(isGene = ifelse(Genes_Slyc != "none", "Gene", "non-Gene")) %>%
  mutate(ismiRNAs = ifelse(miRNAs_Slyc != "none", "miRNA", "non-miRNAs")) %>%
  mutate(isTE = ifelse(TEs_name != "none", "TE", "non-TE")) %>%
  # mutate(isLTRpred = ifelse(TEs_LTRpred != "none", "LTRpred", "non-LTRpredTE")) %>%
  mutate(annotation = "none") %>%
  mutate(annotation = ifelse(isTE == "TE", isTE, annotation)) %>%
  mutate(annotation = ifelse(isGene == "Gene", isGene, annotation)) %>%
  mutate(annotation = ifelse(isPromotor == "Promoter",  "Promoter", annotation)) %>%
  mutate(annotation = ifelse(isGene == "Gene" & isTE == "TE", "TE+Gene", annotation)) %>%
  mutate(annotation = ifelse(isPromotor == "Promoter" & isTE == "TE", "TE+Prom", annotation)) %>%
  mutate(annotation = ifelse(ismiRNAs == "miRNA", ismiRNAs, annotation)) %>%
  mutate(annotation = as.factor(annotation)) %>%
  mutate(alltrue = as.factor(TRUE)) %>%
  dplyr::rename( "#seqname" = seqnames)

# renaming TEs
slydf <- slydf %>%
  mutate(TEs_Order_renamed = fct_recode(TEs_Order,"TIR"="DNA","LINE"="nonLTR", 
                                        "Helitron"="Rolling")) %>%
  mutate(TEs_Order_renamed = fct_collapse(TEs_Order_renamed, "Other_TEs"=c("TE","SSR","Unknown", "PHG"))) %>%
  mutate(TEs_Order_renamed = fct_relevel(TEs_Order_renamed, c("Pararetrovirus", "Helitron", "Other_TEs"), after = Inf))

# adding miRNAs
slydf$miRNAs_Slyc <- sly_bins$miRNAs_Slyc


prefix <- "sRNA_bins200_"
load(file = file.path(  "/projects/TRIBE/srnas/mapped_set3_4_SL30_Penn/", paste0( "sRNA_bins200_", mygenome, "_counts.rdata"))) #mycounts, cnts21, cnts24, 
slydf$srnas_count_21 <- cnts21
slydf$srnas_count_24 <- cnts24
# slydf$srnas_count_all <- cntsall
logcnt2421ratio <- log2(cnts21+1)-log2(cnts24+1)
tmp <- cut(logcnt2421ratio, breaks = c(-Inf,-2,2,Inf))
levels(tmp) <- c("high_24bp", "inbetween", "high_21bp")
slydf$srnas_2124ratio_class <- tmp
tmp <- cut(logcnt2421ratio, breaks = c(-Inf,-1,1,Inf))
levels(tmp) <- c("high_24bp", "inbetween", "high_21bp")
slydf$srnas_2124ratio_relaxed_class <- tmp
table(slydf$srnas_2124ratio_class)
# high_24bp inbetween high_21bp 
#    350145   3700999     89246 
# after correcting meth.diff and adding promoter 
# save(sly_bins, file = file.path(path_comparison, paste0(mygenome, "_bins_plus_statsv7.rdata")))
# after adding sRNA FDR promoter FDR
# save(slydf, file = file.path(path_comparison, paste0(mygenome, "_bins_plus_statsv9.rdata")))

# differnt approach to determine size class, split between 21/22

mymax <- function(x) {y=sum(x); z  <- which(x > y/2); ifelse(length(z) == 1, paste0("high_",z + 19), NA) }

tmp <- slydf %>%
  dplyr::select(matches("srnas_M82C*.*count$")) %>%
  pmap_chr(function(...) mymax(c(...)))
slydf$srnas_M82_size_class  <- tmp

with(slydf, table(srnas_M82_size_class, TEs_Order_renamed))
# EPRVs are mostly 22
with(slydf, table2(srnas_M82_size_class, srnas_2124ratio_class))
with(slydf, table2(srnas_M82_size_class, srnas_dcl2_class2))
                    # srnas_dcl2_class2
# srnas_M82_size_class -D2DESL +D2DESL inbetween non-D2DESL    <NA>
                # 20         9      28     13084      21144       0
                # 21       181      55     27742      43561       0
                # 22      1507      25     40888      54614       0
                # 23        50      58     44797      75666       0
                # 24       921     980    229812     429109       0
                # 25         4       9      4582       7059       0
                # <NA>     461     460    890349    2253235       0
# dcl2 depedent loci are moslty 22

table(tmp+19)

load(file = file.path("/projects/TRIBE/srnas/analysis_bins/sRNA_bins200_SolLyc_bins_raw_processed_dcl2.rdata")) #lrt_dcl2, cnts_dcl2

  tmp2 <- lrt_dcl2$table[, c("logFC", "logCPM", "PValue")] %>%
    mutate(FDR = round(p.adjust(PValue, method = "BH"), 5)) %>%
    mutate(class = cut(FDR, breaks = mybreaks,
                       labels = c("D2DESL", "inbetween", "non-D2DESL"),
                       include.lowest = TRUE)) %>%
    mutate(updown = ifelse(logFC < 0, "-", "+")) %>%
    mutate(class2 = ifelse(class == "D2DESL", paste0(updown, as.character(class)),
                           as.character(class))) %>%
    mutate(logCPM = round(logCPM,2)) %>%
    mutate(logFC = round(logFC,2)) %>%
    dplyr::select(-updown)

slydf$srnas_dcl2_class2 <- tmp2$class2
slydf$srnas_dcl2_class <- tmp2$class
slydf$srnas_dcl2_logFC <- tmp2$logFC

table(slydf$srnas_dcl2_class2)
    # -D2DESL     +D2DESL inbetween  non-D2DESL
    #  3133      1615   1251254   2884388
table(mydf_filter$srnas_dcl2_class2)

  genotype <- sym(paste0("genotype_", myF4))
  mydf_filter <- slydf %>%
    filter(!!genotype == gtfilter)

#------------------------------------ same, loading tmv srna DE data in
  tmp2 <- lrt_tmv$table[, c("logFC", "logCPM", "PValue")] %>%
    mutate(FDR = round(p.adjust(PValue, method = "BH"), 5)) %>%
    mutate(class = cut(FDR, breaks = mybreaks,
                       labels = c("TMV_DESL", "inbetween", "non-D2DESL"),
                       include.lowest = TRUE)) %>%
    mutate(updown = ifelse(logFC < 0, "-", "+")) %>%
    mutate(class2 = ifelse(class == "TMV_DESL", paste0(updown, as.character(class)),
                           as.character(class))) %>%
    mutate(logCPM = round(logCPM,2)) %>%
    mutate(logFC = round(logFC,2)) %>%
    dplyr::select(-updown)

slydf$srnas_tmv_class2 <- tmp2$class2
slydf$srnas_tmv_class <- tmp2$class
slydf$srnas_tmv_logFC <- tmp2$logFC

slydf <- slydf %>%
  mutate(srnas_tmv_class2 = fct_recode(srnas_tmv_class2,"non-TMV_DESL"="non-D2DESL"))

# fix rename issue
colnames(slydf) <- str_replace(colnames(slydf), "_20_class2", "_20bp_counts" )
colnames(slydf) <- str_replace(colnames(slydf), "_21_class2", "_21bp_counts" )
colnames(slydf) <- str_replace(colnames(slydf), "_22_class2", "_22bp_counts" )
colnames(slydf) <- str_replace(colnames(slydf), "_23_class2", "_23bp_counts" )
colnames(slydf) <- str_replace(colnames(slydf), "_24_class2", "_24bp_counts" )
colnames(slydf) <- str_replace(colnames(slydf), "_25_class2", "_25bp_counts" )


with(slydf, table(srnas_tmv_class2, TEs_Order_renamed))
with(slydf, table(srnas_tmv_class2, srnas_dcl2_class2))
with(slydf, table(srnas_tmv_class2, srnas_P4042_class2))
with(slydf, table(srnas_tmv_class2, srnas_P2561_class2))
with(slydf, cor(srnas_tmv_logFC, srnas_dcl2_logFC, method = "spearman"))
# [1] 0.2758596
with(slydf, cor(srnas_tmv_logFC, srnas_P2561_logFC, method = "spearman"))
# [1] 0.05374635

mycounts_dcl2$TEs_Order_renamed <- slydf$TEs_Order_renamed
mycounts_dcl2$annotation <- slydf$annotation
mycounts_dcl2$TEs_Superfamily <- slydf$TEs_Superfamily

colnames(mycounts_dcl2)  <- str_replace(colnames(mycounts_dcl2), "_Sol.*", "")

mycounts_dcl2 %>%
  as_tibble() %>%
  # mutate(across(starts_with("wt"), ~ .x / colSums(.x)))
  mutate(across(contains("-"), ~ .x*100 / sum(.x))) %>%
  # group_by(TEs_Order_renamed) %>%
  group_by(TEs_Superfamily) %>%
  summarize(sum = across(contains("-"), ~ round(sum(.x),1))) %>%
  print(n = Inf)

  # annotation sum$`wt-1` $`wt-2` $`wt-3` $`dclab-1` $`dclab-2` $`wt-1_tmv` $`wt-2_tmv` $`wt-3_tmv`
  # <fct>           <dbl>   <dbl>   <dbl>      <dbl>      <dbl>       <dbl>       <dbl>       <dbl>
# 1 Gene               10      11      10         10         10          12          12          11
# 2 none               13      13      13         13         13          16          17          16
# 3 Promoter           14      14      13         13         16          22          20          20
# 4 TE                 38      39      39         40         40          32          34          34
# 5 TE+Gene             7       7       7          7          7           7           7           7
# 6 TE+Prom            17      16      17         17         15          11          11          11

mycounts_dcl2 %>%
  as_tibble() %>%
  # mutate(across(starts_with("wt"), ~ .x / colSums(.x)))
  mutate(across(contains("-"), ~ .x*100 / sum(.x))) %>%
  # group_by(TEs_Order_renamed) %>%
  group_by(TEs_Order_renamed) %>%
  summarize(sum = across(contains("-"), ~ round(sum(.x),1)))
  # TEs_Order_renamed sum$`wt-1` $`wt-2` $`wt-3` $`dclab-1` $`dclab-2` $`wt-1_tmv` $`wt-2_tmv` $`wt-3_tmv`
  # <fct>                  <dbl>   <dbl>   <dbl>      <dbl>      <dbl>       <dbl>       <dbl>       <dbl>
# 1 TIR                      8.5     7.2     8.4        8.6        6.8         3.9         4.8         5.3
# 2 LTR                     22.8    26.1    24.3       24.3       27.9        26          24.7        24.8
# 3 none                    37.2    38.7    36.6       36.1       38.3        50.3        48.8        47.5
# 4 LINE                     1.2     1       1.2        1.1        1           0.9         1.1         1.1
# 5 Pararetrovirus           1.4     1.5     1.3        0.8        0.6         2.1         2.5         2.3
# 6 Helitron                 0.4     0.7     0.5        0.3        0.9         0.7         0.7         0.7
# 7 Other_TEs               23.6    19.3    22.6       23.6       18.5        10.7        12.3        13.3
# 8 NA                       4.8     5.6     5.2        5          6           5.4         5.1         5.1


# v34 was last to have meth information in
# slydf_meth <- slydf %>%
#     dplyr::select(contains("meth"))

# slydf_PValue <- slydf %>%
#     dplyr::select((ends_with("class") | contains("PValue")))

# v34 was last to have meth information in
# slydf <- slydf %>%
    # dplyr::select(-contains("meth"))
# slydf <- slydf %>%
#     dplyr::select(-(ends_with("class") | contains("PValue")))

save(slydf, file = file.path(path_comparison, paste0(mygenome, "_slydf_plus_statsv38.rdata")))

save(slydf_meth, file = file.path(path_comparison, paste0(mygenome, "_slydf_plus_statsv34_only_meth.rdata")))
save(slydf_PValue, file = file.path(path_comparison, paste0(mygenome, "_slydf_plus_statsv35_only_meth.rdata")))

#}}}
#------------------------------------ dcl2 srnas analysis {{{
with(slydf, table(srnas_dcl2_class2, TEs_Order_renamed))
           # TEs_Order_renamed
# srnas_dcl2_class2     TIR     LTR    none    LINE Pararetrovirus Helitron Other_TEs
  # -DESL         108     506     873      26           1112       19       265
  # +DESL          84     283     908      15              8        3       274
  # inbetween   59144  498545  352898   17956          12252     3212    193718
  # non-DESL   117950 1051313  994203   36959          15081     5704    419755
with(slydf, table2(srnas_dcl2_class2, srnas_P4042_class2))
with(mydf_filter, table2(TEs_Order, srnas_P4042_class2))
with(mydf_filter, table2(TEs_Order, srnas_unique_P4042_class2))
with(mydf_filter, table2(TEs_Order, srnas_dcl2_class2))
with(mydf_filter, table_prop(TEs_Order, srnas_P4042_class2, mymargin = 2))
with(mydf_filter, table_prop(TEs_Order, srnas_P4042_class2, mymargin = 1))
with(slydf, table2(TEs_Order_renamed, srnas_unique_P4042_class2))
with(slydf, table2(srnas_P4042_class2))
with(mydf_filter, table2(srnas_P4042_class2))
with(mydf_filter, prop.table(table(srnas_P4042_class2, deparse.level = 0), margin = 2))

with(mydf_filter, round(prop.table(table(TEs_Order,srnas_P4042_class2, deparse.level = 0), margin = 1)*100,1))
with(mydf_filter, round(prop.table(table(TEs_Order,srnas_unique_P4041_class2, deparse.level = 0), margin = 1)*100,1))

slydf %>%
    dplyr::select(matches("srnas_unique*.*class2$")) %>%
    map(~round(prop.table(table(slydf$TEs_Order,.x, deparse.level = 0), margin = 1)*100,1))

  # examine unique srnas
for (myF4 in F4s_srnas) {
  message(myF4)
  genotype <- sym(paste0("genotype_", myF4))
  mydf_filter <- slydf %>%
    filter(!!genotype == gtfilter)
  tmp <- mydf_filter %>%
    dplyr::select(paste0("srnas_unique_",myF4,"_class2")) %>%
    map(~round(prop.table(table(mydf_filter$TEs_Order,.x, deparse.level = 0), margin = 1)*100,1))
  print(tmp)
}


with(mydf_filter, round(prop.table(table(TEs_Order,srnas_dcl2_class2, deparse.level = 0), margin = 2)*100,1))
with(slydf, round(prop.table(table(TEs_Order_renamed,srnas_dcl2_class2, deparse.level = 0), margin = 2)*100,1))
with(mydf_filter, round(prop.table(table(TEs_Order,srnas_dcl2_class2, deparse.level = 0), margin = 2)*100,1))
with(mydf_filter, round(prop.table(table(TEs_Order,srnas_dcl2_class2, deparse.level = 0), margin = 1)*100,2))

with(mydf_filter, round(prop.table(table(srnas_P4042_class2,srnas_dcl2_class2, deparse.level = 0), margin = 2)*100,1))
with(mydf_filter, round(prop.table(table(srnas_P4042_class2,srnas_dcl2_class2, deparse.level = 0), margin = 1)*100,1))
with(mydf_filter, (table(srnas_P4042_class2)))
with(mydf_filter, (table(srnas_P4042_class2, srnas_dcl2_class2)))

with(slydf, table(TEs_Order_renamed,srnas_dcl2_class2, deparse.level = 0))
                   # -DESL   +DESL inbetween non-DESL
  # TIR                108      84     59144   117950
  # LTR                506     283    498545  1051313
  # none               873     908    352898   994203
  # LINE                26      15     17956    36959
  # Pararetrovirus    1112       8     12252    15081
  # Helitron            19       3      3212     5704
  # Other_TEs          265     274    193718   419755
tmp <- with(slydf, round(prop.table(table(TEs_Order_renamed,srnas_dcl2_class2, deparse.level = 0), margin = 2)*100,2))
                 # -DESL +DESL inbetween non-DESL
  # TIR             3.71  5.33      5.20     4.47
  # LTR            17.39 17.97     43.82    39.81
  # none           30.01 57.65     31.02    37.65
  # LINE            0.89  0.95      1.58     1.40
  # Pararetrovirus 38.23  0.51      1.08     0.57
  # Helitron        0.65  0.19      0.28     0.22
  # Other_TEs       9.11 17.40     17.03    15.89
tmp3 <- with(slydf, round(prop.table(table(TEs_Order_renamed,srnas_dcl2_class2, deparse.level = 0), margin = 1)*100,3))
#                   -DESL  +DESL inbetween non-DESL
#   TIR             0.061  0.047    33.361   66.531
#   LTR             0.033  0.018    32.151   67.798
#   none            0.065  0.067    26.162   73.706
#   LINE            0.047  0.027    32.673   67.252
#   Pararetrovirus  3.908  0.028    43.060   53.003
#   Helitron        0.213  0.034    35.936   63.817
#   Other_TEs       0.043  0.045    31.550   68.363
tmp3 <- with(slydf, round(prop.table(table(TEs_Order_renamed,srnas_dcl2_class2, deparse.level = 0), margin = 1)*100,3))
tmp2 <- as.data.frame(tmp) %>%
  # filter(Var2 != "inbetween") %>%
  # filter(TEs_Order_renamed != "none") %>%
  mutate(TEs_Order_renamed = str_replace(Var1, "none", "non-TE")) %>%
  mutate(TEs_Order_renamed = fct_relevel(TEs_Order_renamed, c("Pararetrovirus", "Helitron", "Other_TEs"), after = Inf)) %>%
  mutate(TEs_Order_renamed = fct_relevel(TEs_Order_renamed, "non-TE"))

gg <- ggplot(tmp2, aes(fill = TEs_Order_renamed, 
                       y = Freq,
                       x = srnas_dcl2_class)) +
# geom_point() +
geom_bar(stat="identity", position="fill") + 
# geom_line() +
theme_bw()+
theme(legend.position = "right",
      axis.text.x = element_text(angle = 300, hjust = 0),
      axis.ticks = element_blank()) +
      scale_fill_manual(values = brewer.pal(7,"Accent"))  +
scale_y_continuous(labels = scales::percent,
                   breaks = scales::pretty_breaks(n = 5))
ggsave("srna_dcl2_dep_vs_TEs_order.pdf",
       gg,
       height = 3,
       width = 3)
#}}}


#------------------------------------ for sara, start here
colnames(slydf)
srna_prefix <- "sRNA_bins200_"
mygenomes <- c("SolLyc", "SolPen")
mygenome <- mygenomes[1]
# load(file = paste0(srna_prefix, mygenome, "_bins_plus_stats.rdata"))
# convert grange into data frame for dplyr magic!!
# load(file = file.path(path_comparison, paste0(mygenome, "sly_bins_plus_statsv20.rdata"))) #slydf
load(file = file.path(path_comparison, paste0(mygenome, "_slydf_plus_statsv31.rdata")))
# load(file = file.path(path_comparison, paste0(mygenome, "_slydf_plus_statsv20.rdata")))

# tmp <- as.data.frame(sly_bins) %>%
#   filter(seqnames %in% c("SL3.0ch11", "SL3.0ch12"))


#  [1] "meth_P1512_CpG_qvalue" "meth_P1512_CHH_qvalue" "meth_P1512_CHG_qvalue" "meth_P1515_CpG_qvalue" "meth_P1515_CHH_qvalue"
#  [6] "meth_P1515_CHG_qvalue" "meth_P2561_CpG_qvalue" "meth_P2561_CHH_qvalue" "meth_P2561_CHG_qvalue" "meth_P3611_CpG_qvalue"
# [11] "meth_P3611_CHH_qvalue" "meth_P3611_CHG_qvalue" "meth_P3612_CpG_qvalue" "meth_P3612_CHH_qvalue" "meth_P3612_CHG_qvalue"
# [16] "meth_P4041_CpG_qvalue" "meth_P4041_CHH_qvalue" "meth_P4041_CHG_qvalue" "meth_P4042_CpG_qvalue" "meth_P4042_CHH_qvalue"
# [21] "meth_P4042_CHG_qvalue"

qvalue_filter <- function(mydata, colname, thres) {
  # mydata <- iris
  # thres <- 5
  # colname <- "Sepal.Length"
  mydata %>%
    filter(!!sym(colname) < thres)
}
# qvalue_filter(iris, "Sepal.Length", 5)

#------------------------------------ saving files

setwd(file.path(path_base, "comparisons_using_big_table/bed_tables/"))
# saving DMRs in bed files:
myvars <- slydf %>%
  dplyr::select(contains("FDR"), contains("qvalue")) %>%
  colnames() 
walk(myvars, ~qvalue_filter(slydf, .x, 0.05) %>% 
     mutate(width = paste0(round(log(!!sym(.x))), TEs_Order)) %>%
     write_tsv( path= paste0("Diffbins200_cutoff_0_05", "_", .x, ".bed"), col_names = T))

# saving list of DEG that are homozygous for each plant
myvars <- slydf %>%
  select(ends_with("FDR")) %>%
  select(contains("rnaseq")) %>%
  colnames() 

# saving bed files of DEGs (e.g. to use with deeptools)
thres <- 0.05
thres2 <- 0.8
for (myF4 in myF4s) {
  # myF4 ="P4042"
  # myF4 ="P2562"
  genotype <- sym(paste0("genotype_", myF4))
  rnaseq <- sym(paste0("rnaseq_", myF4, "_FDR"))

  dgeb <- slydf %>%
    filter(!!rnaseq < thres) %>%
    filter(!!genotype == "00") %>%
    select(Genes_Slyc) %>%
    unique()

  dgebnon <- slydf %>%
    filter(!!rnaseq > thres2) %>%
    filter(!!genotype == "00") %>%
    select(Genes_Slyc) %>%
    unique()

  anno <- genes[["Slyc"]][genes[["Slyc"]]$Name %in% dgeb$Genes_Slyc, ]
  write_tsv(as.data.frame(anno), path= paste0("DEGs_cutoff_0_05_genotype_00", "_", myF4, ".bed"), col_names = F)
  annonon <- genes[["Slyc"]][genes[["Slyc"]]$Name %in% dgebnon$Genes_Slyc, ]
  write_tsv(as.data.frame(annonon), path= paste0("Non_DEGs_cutoff_08_genotype_00", "_", myF4, ".bed"), col_names = F)
}

# saving list of TE bins homozygous for each plant
F4s_all <- c("P1512", "P1515","P2561","P2562","P3611","P3612","P4041","P4042")
gtfilter <- "00"
for (myF4 in F4s_all[-2]) {
  # myF4 <- F4s_all[1]
  genotype <- sym(paste0("genotype_", myF4))
  mydf_filter <- slydf %>%
    filter(!!genotype == gtfilter) %>%
    filter(isTE == "TE")
  write_tsv(as.data.frame(mydf_filter), path= paste0("Transposons_list", "_", myF4, ".csv"), col_names = T)
}

# tmp2 <- slydf[1:1e6,] %>%
tmp2 <- slydf %>%
  group_by(genotype_P4042) %>% 
  # group_modify(~as.data.frame(table(.x$srnas_P4042_class)))
  group_modify(~as.data.frame(100*prop.table(table(.x$srnas_P4042_class))))

for (myF4 in myF4s) {
  for (myF4b in myF4s) {
    # context="CHH"
    if (myF4b == myF4) {
      message("skipping", myF4b)
      next
    }
    # myF4 ="P4041"
    # meth_myplants <-      sym(paste0("meth_", myF4, "_", context, "_qvalue"))
    srna_myplantsclass <- (paste0("srnas_", myF4, "_class"))
    srna_myplantsclass2 <- (paste0("srnas_", myF4b, "_class"))
    genotype <- sym(paste0("genotype_", myF4))
    genotype2 <- sym(paste0("genotype_", myF4b))

    tmp2 <- slydf %>%
      filter(!!genotype == "11") %>%
      filter(!!genotype2 == "11") %>%
      group_by(!!sym(srna_myplantsclass), isPromotor) %>% 
      # group_modify(~as.data.frame(table(.x$srnas_P4042_class)))
      group_modify(~as.data.frame(100*prop.table(table(.x[, srna_myplantsclass2]))))
mytitle <- paste0(srna_myplantsclass2, "_vs_", srna_myplantsclass, ".pdf")
    gg <- ggplot(tmp2, aes(x = !!sym(colnames(tmp2)[1]), 
                           y = Freq, 
                           fill = fct_rev(Var1))) +
# facet_grid(. ~ genotype_P4041 ) +
facet_grid(. ~ isPromotor) +
scale_y_continuous(labels = scales::percent,
                   breaks = scales::pretty_breaks(n = 8)) +
geom_bar(stat="identity",position="fill") +
scale_fill_viridis_d() +
theme(legend.position = "right",
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 300, hjust = 0)) +
 labs(title = mytitle)
ggsave(mytitle,
       gg,
       height = 6,
       width = 6)
  }
}


# for DMR results go DRM_results.R

mybarplot(myF4, myprefix = "srnas_", factor1 = "chromatin_state", factor2 = "isPromotor", gtfilter = "00", mydf = slydf)
mybarplot(myF4, myprefix = "rnaseq_", infix = "", suffix = "_class2", factor1 = "srnas_P4042_class2", factor2 = "TEs_Order", gtfilter = "00", mydf = slydf, mywidth = 30)
mybarplot(myF4, myprefix = "rnaseq_", infix = "", suffix = "_class2", factor1 = "srnas_P4042_class2", factor2 = "chromatin_state", gtfilter = "00", mydf = slydf, mywidth = 30)
mybarplot(myF4, myprefix = "rnaseq_", infix = "", suffix = "_class2", factor1 = "srnas_P4042_class2", factor2 = "alltrue", gtfilter = "00", mydf = slydf, mywidth = 30)
mybarplot(myF4, myprefix = "rnaseq_", infix = "", suffix = "_class2", factor1 = "srnas_P4042_class2", factor2 = "chromatin_state", gtfilter = "00", mydf = slydf, mywidth = 20)

colnames(slydf)
#   [1] "#seqname"                       "start"                          "end"                           
#   [4] "width"                          "strand"                         "ID"                            
#   [7] "genome"                         "type"                           "chromatin_state"               
#  [10] "TEs_name"                       "TEs_Order"                      "TEs_Class"                     
#  [13] "TEs_Superfamily"                "Genes_Slyc"                     "Genes_Spen"                    
#  [16] "Ontology_term_Slyc"             "Note"                           "Promoters_Slyc"                
#  [19] "Promoters_Spen"                 "genotype_M82C"                  "genotype_PenneC"               
#  [22] "genotype_P1512"                 "genotype_P2561"                 "genotype_P2562"                
#  [25] "genotype_P3611"                 "genotype_P3612"                 "genotype_P4041"                
#  [28] "genotype_P4042"                 "srnas_P1512_logFC"              "srnas_P1512_logCPM"            
#  [31] "srnas_P1512_PValue"             "srnas_P1512_FDR"                "srnas_P1512_class"             
#  [34] "srnas_P1512_class2"             "srnas_P2561_logFC"              "srnas_P2561_logCPM"            
#  [37] "srnas_P2561_PValue"             "srnas_P2561_FDR"                "srnas_P2561_class"             
#  [40] "srnas_P2561_class2"             "srnas_P2562_logFC"              "srnas_P2562_logCPM"            
#  [43] "srnas_P2562_PValue"             "srnas_P2562_FDR"                "srnas_P2562_class"             
#  [46] "srnas_P2562_class2"             "srnas_P3611_logFC"              "srnas_P3611_logCPM"            
#  [49] "srnas_P3611_PValue"             "srnas_P3611_FDR"                "srnas_P3611_class"             
#  [52] "srnas_P3611_class2"             "srnas_P3612_logFC"              "srnas_P3612_logCPM"            
#  [55] "srnas_P3612_PValue"             "srnas_P3612_FDR"                "srnas_P3612_class"             
#  [58] "srnas_P3612_class2"             "srnas_P4041_logFC"              "srnas_P4041_logCPM"            
#  [61] "srnas_P4041_PValue"             "srnas_P4041_FDR"                "srnas_P4041_class"             
#  [64] "srnas_P4041_class2"             "srnas_P4042_logFC"              "srnas_P4042_logCPM"            
#  [67] "srnas_P4042_PValue"             "srnas_P4042_FDR"                "srnas_P4042_class"             
#  [70] "srnas_P4042_class2"             "meth_P1512_CpG_qvalue"          "meth_P1512_CpG_meth.diff"      
#  [73] "meth_P1512_CpG_class"           "meth_P1512_CpG_class2"          "meth_P1512_CHH_qvalue"         
#  [76] "meth_P1512_CHH_meth.diff"       "meth_P1512_CHH_class"           "meth_P1512_CHH_class2"         
#  [79] "meth_P1512_CHG_qvalue"          "meth_P1512_CHG_meth.diff"       "meth_P1512_CHG_class"          
#  [82] "meth_P1512_CHG_class2"          "meth_P1515_CpG_qvalue"          "meth_P1515_CpG_meth.diff"      
#  [85] "meth_P1515_CpG_class"           "meth_P1515_CpG_class2"          "meth_P1515_CHH_qvalue"         
#  [88] "meth_P1515_CHH_meth.diff"       "meth_P1515_CHH_class"           "meth_P1515_CHH_class2"         
#  [91] "meth_P1515_CHG_qvalue"          "meth_P1515_CHG_meth.diff"       "meth_P1515_CHG_class"          
#  [94] "meth_P1515_CHG_class2"          "meth_P2561_CpG_qvalue"          "meth_P2561_CpG_meth.diff"      
#  [97] "meth_P2561_CpG_class"           "meth_P2561_CpG_class2"          "meth_P2561_CHH_qvalue"         
# [100] "meth_P2561_CHH_meth.diff"       "meth_P2561_CHH_class"           "meth_P2561_CHH_class2"         
# [103] "meth_P2561_CHG_qvalue"          "meth_P2561_CHG_meth.diff"       "meth_P2561_CHG_class"          
# [106] "meth_P2561_CHG_class2"          "meth_P3611_CpG_qvalue"          "meth_P3611_CpG_meth.diff"      
# [109] "meth_P3611_CpG_class"           "meth_P3611_CpG_class2"          "meth_P3611_CHH_qvalue"         
# [112] "meth_P3611_CHH_meth.diff"       "meth_P3611_CHH_class"           "meth_P3611_CHH_class2"         
# [115] "meth_P3611_CHG_qvalue"          "meth_P3611_CHG_meth.diff"       "meth_P3611_CHG_class"          
# [118] "meth_P3611_CHG_class2"          "meth_P3612_CpG_qvalue"          "meth_P3612_CpG_meth.diff"      
# [121] "meth_P3612_CpG_class"           "meth_P3612_CpG_class2"          "meth_P3612_CHH_qvalue"         
# [124] "meth_P3612_CHH_meth.diff"       "meth_P3612_CHH_class"           "meth_P3612_CHH_class2"         
# [127] "meth_P3612_CHG_qvalue"          "meth_P3612_CHG_meth.diff"       "meth_P3612_CHG_class"          
# [130] "meth_P3612_CHG_class2"          "meth_P4041_CpG_qvalue"          "meth_P4041_CpG_meth.diff"      
# [133] "meth_P4041_CpG_class"           "meth_P4041_CpG_class2"          "meth_P4041_CHH_qvalue"         
# [136] "meth_P4041_CHH_meth.diff"       "meth_P4041_CHH_class"           "meth_P4041_CHH_class2"         
# [139] "meth_P4041_CHG_qvalue"          "meth_P4041_CHG_meth.diff"       "meth_P4041_CHG_class"          
# [142] "meth_P4041_CHG_class2"          "meth_P4042_CpG_qvalue"          "meth_P4042_CpG_meth.diff"      
# [145] "meth_P4042_CpG_class"           "meth_P4042_CpG_class2"          "meth_P4042_CHH_qvalue"         
# [148] "meth_P4042_CHH_meth.diff"       "meth_P4042_CHH_class"           "meth_P4042_CHH_class2"         
# [151] "meth_P4042_CHG_qvalue"          "meth_P4042_CHG_meth.diff"       "meth_P4042_CHG_class"          
# [154] "meth_P4042_CHG_class2"          "rnaseq_P1512_FDR"               "rnaseq_P1512_logFC"            
# [157] "rnaseq_P1512_class2"            "rnaseq_P1512_class"             "rnaseq_P1512_logFC_promoter"   
# [160] "rnaseq_P1512_FDR_promoter"      "rnaseq_P1512_class2_promoter"   "rnaseq_P1512_class_promoter"   
# [163] "rnaseq_P2561_FDR"               "rnaseq_P2561_logFC"             "rnaseq_P2561_class2"           
# [166] "rnaseq_P2561_class"             "rnaseq_P2561_logFC_promoter"    "rnaseq_P2561_FDR_promoter"     
# [169] "rnaseq_P2561_class2_promoter"   "rnaseq_P2561_class_promoter"    "rnaseq_P2562_FDR"              
# [172] "rnaseq_P2562_logFC"             "rnaseq_P2562_class2"            "rnaseq_P2562_class"            
# [175] "rnaseq_P2562_logFC_promoter"    "rnaseq_P2562_FDR_promoter"      "rnaseq_P2562_class2_promoter"  
# [178] "rnaseq_P2562_class_promoter"    "rnaseq_P3611_FDR"               "rnaseq_P3611_logFC"            
# [181] "rnaseq_P3611_class2"            "rnaseq_P3611_class"             "rnaseq_P3611_logFC_promoter"   
# [184] "rnaseq_P3611_FDR_promoter"      "rnaseq_P3611_class2_promoter"   "rnaseq_P3611_class_promoter"   
# [187] "rnaseq_P3612_FDR"               "rnaseq_P3612_logFC"             "rnaseq_P3612_class2"           
# [190] "rnaseq_P3612_class"             "rnaseq_P3612_logFC_promoter"    "rnaseq_P3612_FDR_promoter"     
# [193] "rnaseq_P3612_class2_promoter"   "rnaseq_P3612_class_promoter"    "rnaseq_P4041_FDR"              
# [196] "rnaseq_P4041_logFC"             "rnaseq_P4041_class2"            "rnaseq_P4041_class"            
# [199] "rnaseq_P4041_logFC_promoter"    "rnaseq_P4041_FDR_promoter"      "rnaseq_P4041_class2_promoter"  
# [202] "rnaseq_P4041_class_promoter"    "rnaseq_P4042_FDR"               "rnaseq_P4042_logFC"            
# [205] "rnaseq_P4042_class2"            "rnaseq_P4042_class"             "rnaseq_P4042_logFC_promoter"   
# [208] "rnaseq_P4042_FDR_promoter"      "rnaseq_P4042_class2_promoter"   "rnaseq_P4042_class_promoter"   
# [211] "meth_P1512_CpG_class_stringent" "meth_P1512_CHH_class_stringent" "meth_P1512_CHG_class_stringent"
# [214] "meth_P1515_CpG_class_stringent" "meth_P1515_CHH_class_stringent" "meth_P1515_CHG_class_stringent"
# [217] "meth_P2561_CpG_class_stringent" "meth_P2561_CHH_class_stringent" "meth_P2561_CHG_class_stringent"
# [220] "meth_P3611_CpG_class_stringent" "meth_P3611_CHH_class_stringent" "meth_P3611_CHG_class_stringent"
# [223] "meth_P3612_CpG_class_stringent" "meth_P3612_CHH_class_stringent" "meth_P3612_CHG_class_stringent"
# [226] "meth_P4041_CpG_class_stringent" "meth_P4041_CHH_class_stringent" "meth_P4041_CHG_class_stringent"
# [229] "meth_P4042_CpG_class_stringent" "meth_P4042_CHH_class_stringent" "meth_P4042_CHG_class_stringent"
# [232] "TEs_rnaseq_P1512_FDR"           "TEs_rnaseq_P1512_class2"        "TEs_rnaseq_P1512_logFC"        
# [235] "TEs_rnaseq_P1512_class"         "TEs_rnaseq_P2561_FDR"           "TEs_rnaseq_P2561_logFC"        
# [238] "TEs_rnaseq_P2561_class2"        "TEs_rnaseq_P2561_class"         "TEs_rnaseq_P2562_FDR"          
# [241] "TEs_rnaseq_P2562_logFC"         "TEs_rnaseq_P2562_class2"        "TEs_rnaseq_P2562_class"        
# [244] "TEs_rnaseq_P3611_FDR"           "TEs_rnaseq_P3611_logFC"         "TEs_rnaseq_P3611_class2"       
# [247] "TEs_rnaseq_P3611_class"         "TEs_rnaseq_P3612_FDR"           "TEs_rnaseq_P3612_logFC"        
# [250] "TEs_rnaseq_P3612_class2"        "TEs_rnaseq_P3612_class"         "TEs_rnaseq_P4041_FDR"          
# [253] "TEs_rnaseq_P4041_logFC"         "TEs_rnaseq_P4041_class2"        "TEs_rnaseq_P4041_class"        
# [256] "TEs_rnaseq_P4042_FDR"           "TEs_rnaseq_P4042_logFC"         "TEs_rnaseq_P4042_class2"       
# [259] "TEs_rnaseq_P4042_class"         "isPromotor"                     "isGene"                        
# [262] "isTE"                           "annotation"                     "alltrue"                       
# [265] "srnas_count_21"                 "srnas_count_24"                 "srnas_2124ratio_class"         
# [268] "srnas_2124ratio_relaxed_class" 

mybarplot(myF4s, myprefix = "srnas_", overlaptype = "chromatin_state", gtfilter = "00")
mybarplot(myF4s, myprefix = "meth_", infix = "_CHG", overlaptype = "chromatin_state", gtfilter = "00")
mybarplot(myF4s, myprefix = "meth_", infix = "_CHH", overlaptype = "chromatin_state", gtfilter = "00")
mybarplot(myF4s, myprefix = "meth_", infix = "_CpG", overlaptype = "chromatin_state", gtfilter = "00")
mybarplot(myF4s, myprefix = "meth_", infix = "_CHG", overlaptype = "TEs_Order", gtfilter = "00")
mybarplot(myF4s, myprefix = "meth_", infix = "_CHH", overlaptype = "TEs_Order", gtfilter = "00")
mybarplot(myF4s, myprefix = "meth_", infix = "_CpG", overlaptype = "TEs_Order", gtfilter = "00")
mybarplot(myF4s, myprefix = "", infix = "", overlaptype = "TEs_Order", gtfilter = "00", myplantsclass = "alltrue")


tmp3 <- tmp %>%
  filter(Genes_Slyc %in% c("Solyc08g028970.3.1"))

mybreaks <- c(0,0.1,0.9,1)
myF4s <- c("P1512", "P2561","P3611","P3612","P4041", "P4042")
for (myF4 in myF4s) {
  for (context in c("CpG", "CHH", "CHG")) {
    # context="CHH"
    # myF4 ="P4041"
    meth_myplants <-      sym(paste0("meth_", myF4, "_", context, "_qvalue"))
    meth_myplantsclass <- sym(paste0("meth_", myF4, "_", context, "_class"))
    tmp2 <- tmp2 %>%
      dplyr::mutate(!!meth_myplantsclass := cut(!!meth_myplants, breaks = mybreaks))
  }
}

mytab <- function(tab) round(prop.table(table(tab)),2)

slydf %>%
  select(ends_with("class")) %>%
  map(mytab)


tmp2 %>%
  select(ends_with("class"), chromatin_state) %>%
  map(mytab)

  # tmp <- tmp2 %>%
  #   select(!!meth_myplantsclass) %>%
  #   table()



table((tmp2$meth_P1512_CpG_qvalue))
################loop to do a plot per plant!!!!!!###thanks to Seb. 
myF4s <- c("P1512", "P2561", "P2562","P3611","P3612","P4041", "P4042")
for (myF4 in myF4s) {
  # myF4 = myF4s[1]
  logFC_srnas_myplants <- sym(paste0("srnas_", myF4, "_logFC"))
  genotype_myplants <- sym(paste0("genotype_", myF4))
  
  message(logFC_srnas_myplants)
  gg <-
    ggplot(tmp, aes(
      x = start,
      y = !!logFC_srnas_myplants,
      color = !!genotype_myplants
      )) +
    geom_point() +
    facet_grid(seqnames ~ .) +
    theme(legend.title = element_blank(),
          strip.text.y = element_text(size = 12, angle = 0)) +
    guides(fill = guide_legend(title = NULL)) +
    scale_colour_discrete(labels = c("Slyc homozygous", "heterozygous", "Spenn homozygous"))
  ggsave(paste0("sRNA_cov_", myF4, ".png"),
         gg,
         height = 15,
         width = 49)
}

tmp2 <- tmp %>%
  filter(Genes_Slyc != "none")

gg <- ggplot(tmp, aes(x = start, y = srnas_P4042_logCPM, color = genotype_10.P4042.1_MapppedOn_SL30)) +
  geom_line() +
  facet_grid(seqnames ~ .)

gg <- ggplot(tmp, aes(x = start, y = as.numeric(meth_P4042_CpG_meth.diff), color = genotype_P4042)) +
  geom_line() +
  facet_grid(seqnames ~ .)
ggsave(file=file.path("methP4042_diff_cpg.pdf"),gg,height=10,width=30)

gg <- ggplot(tmp, aes(x = start, y = as.numeric(meth_P4042_CHG_meth.diff), color = genotype_P4042)) +
  geom_line() +
  facet_grid(seqnames ~ .)
ggsave(file=file.path("methP4042_diff_CHG.pdf"),gg,height=10,width=30)

gg <- ggplot(tmp, aes(x = start, 
                      y = as.numeric(srnas_P4042_logFC), 
                      color = genotype_10.P4042.1_MapppedOn_SL30)) +
  geom_point() +
  ylim(-5, 5) +
  facet_grid(seqnames ~ .)

ggsave(file=file.path(path_srna,"logFC.pdf"),gg,height=5,width=30)





p11<-ggplot(tmp, aes(x=-(!!meth_myplants_CHH), y=!!logFC_srnas_myplants)) +
  geom_vline(xintercept = 0,colour= "grey")+
  geom_hline(yintercept = 0, colour= "grey")+
  theme_bw()+
  geom_point(shape=19, aes(color = TEs_Superfamily))
  # scale_color_viridis(option = "magma")

# library(magg
tmp %$% 
  cor(meth_P3611_CHH_meth.diff, srnas_P3611_logFC , use = "complete.obs")

cor((!!meth_myplants_CHH), !!( logFC_srnas_myplants ))

#p11
library(viridis)
options(bitmapType='cairo')
ggsave(paste0("all_methCHH_", myF4, "_sRNAs_scatterplot.png"), p11, height=30,width=30)

## troubleshooting
 filter(tmp, str_detect(Genes_Slyc, "Solyc10g054480"))
