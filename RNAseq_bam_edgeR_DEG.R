library(readr)
library("Rsamtools")
library("Rsubread")
library(edgeR)
library(dplyr)
library(stringr)
library(purrr)
library(rtracklayer)
source("/projects/TRIBE/scripts/project_settings.R")

# load("/data/public_data/tomato/additional_resources/R_code_and_objects/R-objects-bins-SL30-Penn_v6_annot.rdata")
load(file = file.path(path_comparison, paste0("R-objects-bins-SL30-Penn_v6_annot.rdata")))
(myPaths)
#                                                                                                                    genome 
#                                          "/data/public_data/tomato/merged_genomes/genomes_SL30_penn_organelles_merged.fa" 
#                                                                                                                genes_Slyc 
#                                                     "/data/public_data/tomato/ITAG3.2/ITAG3.2_gene_models.gff_sorted.gff" 
#                                                                                                                genes_Spen 
# "/data/public_data/tomato/Solanum_pennellii_complete/Solanum_pennellii_annotation/spenn_v20_gene_models_annot_sorted.gff" 
#                                                                                                                   TEsSlyc 
#                                                    "/data/public_data/tomato/ITAG3.2/ITAG3.2_REPET_repeats_agressive.gff" 

# DEG gene analysis from BAM files using edgeR
# /projects/TRIBE/RNA_seq/20170901_RNAseq_parents_and_F4_preprocessed_by_Sara/mapped_data/STAR_SL3.0_plus_organelles
# /projects/TRIBE/RNA_seq/20200217_RNAseq_parents_and_F4_seb_analyis_sara_mapped
prefix <- "RNAseq_TEs_"
prefix <- "RNAseq_genes_"

mybreaks <- c(0, 0.05, 0.9, 1)
path_base_rnaseq <-  "/projects/TRIBE/RNA_seq/"
path_map  <- file.path(path_base_rnaseq, "20170901_RNAseq_parents_and_F4_preprocessed_by_Sara/mapped_data/STAR_SL3.0_plus_organelles")
mypath <- "/projects/TRIBE/srnas/mapped_set3_4_SL30_Penn"
F4s_all <- c("P1512", "P1515","P2561","P2562","P3611","P3612","P4041","P4042")
F4s_rnaseq <- F4s_all[! F4s_all %in% "P1515"]
genes_sly <- "/data/public_data/tomato/ITAG3.2/ITAG3.2_gene_models.gtf2"
# selgen <- read.delim("/projects/TRIBE/RNA_seq/10selectedgenes_data.csv")
TEs_path <- "/data/public_data/tomato/ITAG3.2/ITAG3.2_REPET_repeats_agressive.gtf"

fls <- list.files(path_map, pattern="out.bam$", full.names =F)
# fls_full <- list.files(path_map, pattern="out.bam$", full.names =T)

  meta <- data.frame(File = fls) %>%
    filter(!str_detect(File, "DDM")) %>%
    mutate(id = str_split(File,"-") %>% map_chr(.,1)) %>%
    mutate(lib = str_split(File,"-") %>% map_chr(.,2)) %>%
    mutate(rep = str_split(File,"-") %>% map_chr(.,3) %>% substr(1,1))

setwd(file.path(path_base_rnaseq, "20200217_RNAseq_parents_and_F4_seb_analyis_sara_mapped"))
write_csv(meta, "RNAseq_metadata.csv")

# using Rsubread to import read-count matrix
fc <- featureCounts(files=file.path(path_map, meta$File), 
                    annot.ext=genes_sly, 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE)

fcTEs <- featureCounts(files=file.path(path_map, meta$File), 
                    annot.ext=TEs_path, 
                    GTF.attrType="transcript_id",
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE)

save(fc, fcTEs, file = file.path(path_base_rnaseq, paste0(prefix, "SolLyc_edgeR_subread_objects.rdata")))

myPaths["genes_Slyc"]
rna_cnts <- list()
tes_cnts <- list()

myparent <- "PenneC"
# myparent <- "M82C"
for (myF4 in F4s_rnaseq) {
  message(myF4)
  # myF4 <-  myF4s[7]
  # select only parents and F4s for time being:

  meta_sub <- meta %>%
    dplyr::mutate(sub = lib == myparent | lib == myF4)
  group <- subset(meta_sub, sub)$lib
  # counts_sub <- fc$counts[,which(meta_sub$sub)]
  counts_sub <- fcTEs$counts[,which(meta_sub$sub)]
  dgList <- DGEList(counts=counts_sub,
                    group=group,
                    remove.zeros = F,
                    genes=rownames(fcTEs$counts)
                    # genes=rownames(fc$counts)
  )
  dgList <- calcNormFactors(dgList, method="TMM")
  design <- model.matrix(~group)
  dgList <- estimateDisp(dgList, design)
  #   png(file=paste0(prefix, myF4, "_BCV", ".png"), width=1000, height=1000, type = "cairo-png")
  #   plotBCV(dgList)
  #   dev.off()
  fit <- glmFit(dgList,design)
  lrt <- glmLRT(fit, coef=2)

  lrt$table <- lrt$table %>%
    mutate(FDR = round(p.adjust(PValue, method = "BH"), 5)) %>%
    mutate(class = cut(FDR, breaks = mybreaks,
                       labels = c("DEG", "inbetween", "non-DEG"),
                       include.lowest = TRUE)) %>%
    mutate(updown = ifelse(logFC < 0, "-", "+")) %>%
    mutate(class2 = ifelse(class == "DEG", paste0(updown, as.character(class)),
                           as.character(class))) %>%
    mutate(logCPM = round(logCPM,2)) %>%
    mutate(logFC = round(logFC,2))

  # png(file=paste0(prefix, myF4, "DE", ".png"), width=1000, height=1000, type = "cairo-png")
  # plotSmear(lrt, de.tags=deGenes)
  # dev.off()
  # tags <- topTags(lrt, n = 1e5, p.value = 0.05)
  table_all <- cbind(lrt$table, counts_sub, fcTEs$annotation)
  # table_all <- cbind(lrt$table, counts_sub, fc$annotation)

  table_all <- table_all %>%
    # mutate(genename = str_replace(GeneID, "gene:", "")) %>%
    # mutate(genename = str_replace(genename, "[.][0-9]", "")) %>%
    mutate(start = str_split(as.character(Start),";") %>% map_chr(., 1)) %>%
    mutate(end = str_split(as.character(End),";") %>% map_chr(., dplyr::last)) %>%
    mutate(seqnames = str_split(as.character(Chr),";") %>% map_chr(., 1)) %>%
    mutate(strand = str_split(as.character(Strand),";") %>% map_chr(., 1)) %>%
    # mutate(strand = substr(as.character(X3), 10, 10)) %>%
    mutate(FDR = p.adjust(PValue, method="BH")) %>%
    dplyr::select(-c(Start, Strand, End, Chr))

  # for TEs... uncomment
   id <- table_all$GeneID
   tmp2 <- anno[(id),]
   table_all <- cbind(table_all, tmp2)

  lrt$table_all <- table_all
  rna_cnts[[myF4]] <- lrt
  write_csv(table_all, path = paste0(prefix, "SolLyc", "_", myF4, "_vs_", myparent, "_table.csv"))
}

# saveRDS(rna_cnts, file = file.path(path_base_rnaseq, paste0(prefix, "SolLyc_edgeR_Robject.rds")))
saveRDS(rna_cnts, file = file.path(path_base_rnaseq, paste0(prefix, "_vs_", myparent,  "_SolLyc_edgeR_Robject.rds")))
# P1512
#      -DEG      +DEG inbetween   non-DEG
#      1678      1726     14994     17370
# P2561
#      -DEG      +DEG inbetween   non-DEG
#       453       612     12386     22317
# P2562
#      -DEG      +DEG inbetween   non-DEG
#      2911      2698     16734     13425
# P3611
#      -DEG      +DEG inbetween   non-DEG
#      1590      1338     15734     17106
# P3612
#      -DEG      +DEG inbetween   non-DEG
#      1422      1124     15709     17513
# P4041
#      -DEG      +DEG inbetween   non-DEG
#      1057       980     15072     18659
# P4042
#      -DEG      +DEG inbetween   non-DEG
#       740       815     13340     20873
# Penn:
map(rna_cnts, ~ table(.x$table_all$class2))

anno=as.data.frame(TEs[[1]])[,-c(1:9)]
rownames(anno)=anno$ID 
table(lrt$table_all$class)

#------------------------------------ dcl2 RNA-Seq analysis
path_map_dcl2 <- file.path(path_base, "RNA_seq/dcl2_rnaseq/mapped")
meta_dcl2 <- read_csv(file.path(path_base, "RNA_seq/dcl2_rnaseq/RNA_seq_metadata_dcl2.csv"))

# only select non infected libs
meta_dcl2_nonvirus <- meta_dcl2  %>%
  filter(virus_infected == "NO")

fc <- featureCounts(files=file.path(path_map_dcl2, meta_dcl2_nonvirus$File), 
                    annot.ext=genes_sly, 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE)

  group <- meta_dcl2_nonvirus$plant
  # counts_sub <- fc$counts[,which(meta_sub$sub)]
  dgList <- DGEList(counts=fc$counts,
                    group=group,
                    remove.zeros = F,
                    genes=rownames(fc$counts)
  )
  dgList <- calcNormFactors(dgList, method="TMM")
  design <- model.matrix(~rev(group))
  dgList <- estimateDisp(dgList, design)
  #   png(file=paste0(prefix, myF4, "_BCV", ".png"), width=1000, height=1000, type = "cairo-png")
  #   plotBCV(dgList)
  #   dev.off()
  fit <- glmFit(dgList,design)
lrt_dcl2 <- glmLRT(fit, coef=2)


lrt_dcl2$table <- lrt_dcl2$table %>%
  mutate(FDR = round(p.adjust(PValue, method = "BH"), 5)) %>%
  mutate(class = cut(FDR, breaks = mybreaks,
                       labels = c("DEG", "inbetween", "non-DEG"),
                       include.lowest = TRUE)) %>%
  # mutate(class_stringent2 = ifelse(abs(logFC) > 2 && class_stringent == "DEG", "DEG", "class_stringent")) %>%
  mutate(class = ifelse(abs(logFC) < 2 & class == "DEG", "inbetween", as.character(class))) %>%
  mutate(updown = ifelse(logFC < 0, "-", "+")) %>%
  mutate(class2 = ifelse(class == "DEG", paste0(updown, as.character(class)),
                         as.character(class))) %>%
  mutate(logCPM = round(logCPM,2)) %>%
  mutate(logFC = round(logFC,2))

  table_all <- cbind(lrt_dcl2$table, fc$counts, fc$annotation)

lrt_dcl2$table_all <- table_all %>%
    mutate(start = str_split(as.character(Start),";") %>% map_chr(., 1)) %>%
    mutate(end = str_split(as.character(End),";") %>% map_chr(., dplyr::last)) %>%
    mutate(seqnames = str_split(as.character(Chr),";") %>% map_chr(., 1)) %>%
    mutate(strand = str_split(as.character(Strand),";") %>% map_chr(., 1)) %>%
    # mutate(strand = substr(as.character(X3), 10, 10)) %>%
    mutate(FDR = p.adjust(PValue, method="BH")) %>%
    dplyr::select(-c(Start, Strand, End, Chr)) %>%
  mutate_at(vars(matches("*class*")),
           ~fct_recode(., "-D2G" = "-DEG",
                          "+D2G" = "+DEG",
                          "non-D2G" = "non-DEG",
                          "D2G" = "DEG" ))

prefix <- "RNAseq_dcl2b_nonvirus"
write_csv(table_all, path = paste0(prefix, "SolLyc", "_table.csv"))

saveRDS(lrt_dcl2, file = file.path(path_base_rnaseq, paste0(prefix, "SolLyc_edgeR_Robject_dcl2_deg.rds")))

deGenes <- rownames(lrt_dcl2)[lrt_dcl2$table$class == "DEG"]
png(file=paste0(prefix, "dcl2", "_BCV", ".png"), width=1000, height=1000, type = "cairo-png")
plotBCV(dgList)
dev.off()
png(file=paste0(prefix, "dcl2", "DE", ".png"), width=1000, height=1000, type = "cairo-png")
plotSmear(lrt_dcl2, de.tags=deGenes)
dev.off()
table(lrt_dcl2$table$class2)

      # -DEG -inbetween       +DEG +inbetween  inbetween    non-DEG
      # 1294       6895        732       6859       7544      12444
table(lrt_dcl2$table$class)
#       DEG inbetween   non-DEG 
#     15780      7544     12444 
table(lrt_dcl2$table$class2)
#      -DEG      +DEG inbetween   non-DEG 
#      8189      7591      7544     12444 
table(lrt_dcl2$table$class_stringent)
#       DEG inbetween   non-DEG 
#     12618     10706     12444 

## non-virus only
table(lrt_dcl2$table$class)
      # DEG inbetween   non-DEG
     # 1777     20554     13437
table(lrt_dcl2$table$class2)
     # -DEG      +DEG inbetween   non-DEG
     #  895       882     20554     13437


# supp: It will be good to show that DCL2 are not DE (or not so much) in our F4s
dcl2genes <- c("Solyc06g048960.3.1","Solyc07g005030.3.1","Solyc08g067210.3.1","Solyc10g005130.3.1","Solyc11g008520.2.1","Solyc11g008530.2.1","Solyc11g008540.2.1")

tmp <- slydf %>%
  filter(Genes_Slyc %in% dcl2genes) %>%
  distinct(Genes_Slyc, .keep_all= TRUE) %>%
  dplyr::select(c(1:28, contains("rnaseq"), c("srnas_dcl2_class2"))) %>%
  dplyr::select(!(contains("meth") | ends_with("count") | starts_with("TEs") | contains("promoter") |contains("unique") | contains("srnas") | contains("class3") | ends_with("class") | contains("class_")))

write_csv(tmp, "rnaseq_dcl2genes_F4s_dcl2lib_bigtable.csv")
