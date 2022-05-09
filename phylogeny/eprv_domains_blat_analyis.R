path_base <-  "/projects/TRIBE"
path_giri <- "/projects/TRIBE/EPRV_annotation/EPRV_giri"
setwd(path_giri)
source(file.path(path_base, "scripts/project_settings.R"))
source(file.path(path_base, "scripts/functions_phasing.R"))
source("~/code/R-code-misc-seb/R-functions-seb.r")
library("hiReadsProcessor")

mypsl <- read.psl("ORF_ntseq_fromGIRI_noheader.psl", bestScoring = F, asGRanges = T, removeFile = F)
eprv_domains_blat <- read.psl("ORF_ntseq_fromGIRI_relaxed.nohead.psl", bestScoring = F, asGRanges = T, removeFile = F)
eprv_domains_blat_penn <- read.psl("ORF_ntseq_fromGIRI_penn_relaxed.nohead.psl", bestScoring = F, asGRanges = T, removeFile = F)
eprv_blat_penn <- read.psl("LycEPRV_I_plus_subsections.fa_penn_relaxed.nohead.psl", bestScoring = F, asGRanges = T, removeFile = F)
eprv_blat <- read.psl("LycEPRV_I_plus_subsections.fa_relaxed.nohead.psl", bestScoring = F, asGRanges = T, removeFile = F)

# containing full blat EPRV results!!:
eprv_blat_full_penn <- eprv_blat_penn[eprv_blat_penn@elementMetadata$qName == "LycEPRV_I",]
eprv_blat_full <- eprv_blat[eprv_blat@elementMetadata$qName == "LycEPRV_I",]

  penn <- eprv_domains_blat_penn[eprv_domains_blat_penn@elementMetadata$qName == mytype,]

TEs_path <- "/data/public_data/tomato/ITAG3.2/ITAG3.2_REPET_repeats_agressive.gtf"
load("/data/public_data/tomato/additional_resources/R_code_and_objects/R-objects-bins-SL30-Penn_v6_annot.rdata")
# load(file = file.path(path_comparison, paste0(mygenome, "_slydf_plus_statsv30.rdata")))
load(file = file.path(path_comparison, paste0(mygenome, "_bins_plus_statsv23.rdata")))
head(slydf)
TEs <- TEs[["Slyc"]]
sly_bins$srnas_dcl2_class2 <- slydf$srnas_dcl2_class2
sly_bins$srnas_tmv_class2 <- slydf$srnas_tmv_class2

    tmp <- overlap_annot(mypsl, TEs, "anno", "Class")
    mypsl <- overlap_annot(mypsl, TEs, "anno", "Order")
    mypsl <- overlap_annot(mypsl, TEs, "anno2", "Superfamily")
    eprv_domains_blat <- overlap_annot(eprv_domains_blat, TEs, "Order", "Order")
    eprv_domains_blat <- overlap_annot(eprv_domains_blat, TEs, "Superfamily", "Superfamily")
    eprv_domains_blat <- overlap_annot(eprv_domains_blat, sly_bins, "srnas_P1512_class2", "srnas_P1512_class2")
    eprv_domains_blat <- overlap_annot(eprv_domains_blat, sly_bins, "srnas_P4042_class2", "srnas_P4042_class2")
    eprv_domains_blat <- overlap_annot(eprv_domains_blat, sly_bins, "srnas_P4042_class2", "srnas_P4042_class2")
    eprv_domains_blat <- overlap_annot(eprv_domains_blat, sly_bins, "srnas_tmv_class2", "srnas_tmv_class2")

    eprv_blat_full <- overlap_annot(eprv_blat_full, sly_bins, "srnas_P4042_class2", "srnas_P4042_class2")
    eprv_blat_full <- overlap_annot(eprv_blat_full, sly_bins, "TEs_Order_renamed", "TEs_Order_renamed")
    eprv_blat_full <- overlap_annot(eprv_blat_full, sly_bins, "srnas_P1512_class2", "srnas_P1512_class2")
    eprv_blat_full <- overlap_annot(eprv_blat_full, sly_bins, "srnas_tmv_class2", "srnas_tmv_class2")
    eprv_blat_full <- overlap_annot(eprv_blat_full, sly_bins, "srnas_dcl2_class2", "srnas_dcl2_class2")

with(as.data.frame(eprv_domains_blat), table(Superfamily))
with(as.data.frame(eprv_domains_blat), table(srnas_P1512_class2))
with(as.data.frame(eprv_domains_blat), table(qName))
with(as.data.frame(eprv_domains_blat), table2d(Superfamily, qName))
with(as.data.frame(eprv_domains_blat), table2d(Superfamily, srnas_P1512_class2))
with(as.data.frame(eprv_domains_blat), table2d(qName, srnas_P1512_class2))
with(as.data.frame(eprv_domains_blat), table2d(anno2, qName))
             # qName  env inclusion_body movement_protein  pol   NA  Sum
# anno2
# Confused_TE           1              0                0    0    0    1
# Copia                 1              0                1    0    0    2
# EPRV                347            207              364  375    0 1293
# Gypsy                 3              1                4    0    0    8
# none                  0              0                1    0    0    1
# SSR                   1              0                1    0    0    2
# Unclassified          1              0                0    0    0    1
# NA                   36             13               37   31    0  117
# Sum                 390            221              408  406    0 1425

with(as.data.frame(mypsl), table(anno))
with(as.data.frame(mypsl), table(qName))
with(as.data.frame(mypsl), table2d(anno, qName))
with(as.data.frame(mypsl), table(anno, qName))
                # qName
# anno             env inclusion_body movement_protein pol
  # LTR              2              0                1   1
  # none            54             27               82 102
  # Pararetrovirus 120             71               71 141
  # TE               1              0                0   0
save(eprv_domains_blat, file = file.path(path_giri, "eprv_domains_blat.rdata"))

#------------------------------------ pylogeny
# Extract fasta of all blat findings and do multiple alignment
# Name of fasta entries should reflect their annotation, i.e. overlap with Superfamily, DESL etc.
names(eprv_domains_blat) <- with(as.data.frame(eprv_domains_blat), 
                                 paste(qName, seqnames, start, Superfamily, srnas_P1512_class2, collpase = "", sep = "_"))

names(eprv_domains_blat) <- with(as.data.frame(eprv_domains_blat), 
                                 paste(qName, seqnames, start, Superfamily, "P4042", srnas_P4042_class2, "D2SL",srnas_dcl2_class2, "TMV", srnas_tmv_class2, collpase = "", sep = "_")) # 

names(eprv_blat_full) <- with(as.data.frame(eprv_blat_full), 
                                 paste(qName, seqnames, start, "P4042", srnas_P4042_class2, "D2SL",srnas_dcl2_class2, "TMV", srnas_tmv_class2, collpase = "", sep = "_")) # 

names(eprv_blat_full_penn) <- with(as.data.frame(eprv_blat_full_penn), 
                                 paste(seqnames, start, "Penn", collpase = "", sep = "_")) # 

# extract sequences
tomato_fa <- readDNAStringSet(genome_sly_path, format="fasta")
names(tomato_fa) <-  str_replace(names(tomato_fa)," ","")

penn_fa <- readDNAStringSet(genome_penn_path, format="fasta")
# [1] "ENA|HG975451|HG975451.1" ...
names(penn_fa) <-  str_replace(names(penn_fa)," .*","")

types <- unique(eprv_domains_blat@elementMetadata$qName)
table(eprv_domains_blat@elementMetadata$qName)
             # env   inclusion_body movement_protein              pol
             # 390              221              408              406
table(eprv_domains_blat_penn@elementMetadata$qName)
             # env   inclusion_body movement_protein              pol
             # 361              141              360              360

with(as.data.frame(eprv_domains_blat), table(srnas_dcl2_class2, qName))
with(as.data.frame(eprv_domains_blat), table(srnas_tmv_class2, qName))
with(as.data.frame(eprv_domains_blat), table(srnas_tmv_class2, srnas_dcl2_class2))
with(slydf, table(srnas_tmv_class2, srnas_dcl2_class2, TEs_Order_renamed))

, , TEs_Order_renamed = TIR

                srnas_dcl2_class2
srnas_tmv_class2 -D2DESL +D2DESL inbetween non-D2DESL
      -TMV_DESL       70       3      2918       5020
      +TMV_DESL        5      33      1164       1604

, , TEs_Order_renamed = Pararetrovirus

                srnas_dcl2_class2
srnas_tmv_class2 -D2DESL +D2DESL inbetween non-D2DESL
      -TMV_DESL       25       1        65         97
      +TMV_DESL      761       1      3012       1376


for (mytype in types) {
  # mytype <- "pol"
  pol <- eprv_domains_blat[eprv_domains_blat@elementMetadata$qName == mytype,]
  penn <- eprv_domains_blat_penn[eprv_domains_blat_penn@elementMetadata$qName == mytype,]
  tmp <- tomato_fa[pol,]
  seqs_penn <- penn_fa[penn,]
  # writeXStringSet(tmp, paste0(mytype, "_P4042_blat_extracted_seqs.fa"), append=FALSE,
  #                 compress=FALSE, compression_level=NA, format="fasta")
  # penn
  writeXStringSet(seqs_penn, paste0(mytype, "_penn_blat_extracted_seqs.fa"), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")
}

tomato_fa <- readDNAStringSet(genome_sly_path, format="fasta")
names(tomato_fa) <-  str_replace(names(tomato_fa)," ","")
# extract sequences

# full EPRV
mytype <- "full_eprvs"

  tmp <- tomato_fa[eprv_blat_full,]
  seqs_penn <- penn_fa[eprv_blat_full_penn,]
  writeXStringSet(tmp, paste0(mytype, "_lyc_blat_extracted_seqs.fa"), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")
  writeXStringSet(seqs_penn, paste0(mytype, "_penn_blat_extracted_seqs.fa"), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")


eprv_df <- as.data.frame(eprv_blat_full)
eprv_df$size <- with(eprv_df, qEnd - qStart)

eprv_df$size_class <- cut(eprv_df$size, 4)
with(eprv_df, table2(srnas_P4042_class2, size_class))
with(eprv_df, table2(srnas_dcl2_class2, size_class))
with(eprv_df, table2(srnas_tmv_class2, size_class))
with(eprv_df, table2(srnas_P1512_class2, size_class))

srnas_dcl2_class2
# run t_coffee on the results! bash:

t_coffee movement_protein_blat_extracted_seqs.fa

file=movement_protein_blat_extracted_seqs.fa
file=inclusion_body_blat_extracted_seqs
trimal -in ${file}.aln -out ${file}_gappy.aln -gappyout

