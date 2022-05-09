path_base <-  "/projects/TRIBE"
path_giri <- "/projects/TRIBE/EPRV_annotation/EPRV_giri"
source(file.path(path_base, "scripts/project_settings.R"))
source(file.path(path_base, "scripts/functions_phasing.R"))
source("~/code/R-code-misc-seb/R-functions-seb.r")
library("hiReadsProcessor")

mypsl <- read.psl("ORF_ntseq_fromGIRI_noheader.psl", bestScoring = F, asGRanges = T, removeFile = F)
eprv_domains_blat <- read.psl("ORF_ntseq_fromGIRI_relaxed.nohead.psl", bestScoring = F, asGRanges = T, removeFile = F)

TEs_path <- "/data/public_data/tomato/ITAG3.2/ITAG3.2_REPET_repeats_agressive.gtf"
load("/data/public_data/tomato/additional_resources/R_code_and_objects/R-objects-bins-SL30-Penn_v6_annot.rdata")
TEs <- TEs[["Slyc"]]

    tmp <- overlap_annot(mypsl, TEs, "anno", "Class")
    mypsl <- overlap_annot(mypsl, TEs, "anno", "Order")
    mypsl <- overlap_annot(mypsl, TEs, "anno2", "Superfamily")
    eprv_domains_blat <- overlap_annot(eprv_domains_blat, TEs, "anno", "Order")
    eprv_domains_blat <- overlap_annot(eprv_domains_blat, TEs, "anno2", "Superfamily")

with(as.data.frame(eprv_domains_blat), table(anno))
with(as.data.frame(eprv_domains_blat), table(qName))
with(as.data.frame(eprv_domains_blat), table2d(anno, qName))
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

