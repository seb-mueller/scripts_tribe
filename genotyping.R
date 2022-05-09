
## bash preprocessing of 
mummer=/applications/MUMmer/MUMmer3.23/
#delta filter
for file in $(ls *[^r].delta); do
        echo ${file%.delta}.filter.delta
        if [ ! -e "${file%.delta}.filter.delta" ]
        then
    $mummer/delta-filter -r -q $file > ${file%.delta}.filter.delta &
        fi
done

#creating coord, snps and its counterparts gff and vcf:
for file in $(ls Solanum*[^r].filter.delta); do
# for file in $(ls Solanum*M82*.filter.delta); do
        echo $file
        $mummer/show-coords -rcl $file > $file.coords
        cat $file.coords | tail -n +6 | perl -lane '{print "$F[17]\tNucmer_delta_filter\thomolog\t$F[0]\t$F[1]\t\.\t\.\t\.\tID=@F"}' > $file.coords.gff3
        $mummer/show-snps -CTr $file > $file.snps
        cat /home_old/sm934/analysis/tomato/mummer/header.vcf > ${file}.vcf
        cat ${file}.snps | tail -n +5 | perl -lane 'if ($F[1]=~/[ACGT]/ & $F[2]=~/[ACGT]/) {print "$F[8]\t$F[0]\t\.\t$F[1]\t$F[2]\t10\tPASS\tAF=1"}' >> ${file}.vcf
        /applications/IGVTools/v2_3_53/igvtools index ${file}.vcf
done



## using vcftools:

# Output allele frequency for all sites in the input vcf file from chromosome 1
# /data/public_data/tomato/assembly_build_3.00/inhouse_annotation/Solanum_Heinz30_Penn.filter.delta.vcf
snps=Solanum_Heinz30_Penn.filter.delta.vcf
ref=/data/public_data/tomato/S_lycopersicum/S_lycopersicum_chromosomes.3.00.fa

cd /projects/TRIBE/RNA_seq/20170901_RNAseq_parents_and_F4_preprocessed_by_Sara/mapped_data/STAR_SL3.0_plus_organelles
bam=22-P4042-2_MapppedOn_SL30_Aligned.sortedByCoord.out.bam
bam=9-P4041-1_MapppedOn_SL30_Aligned.sortedByCoord.out.R.vcf
for bam in *.out.bam; do
  echo "Processing $bam file..."
  bam=${bam%.bam}
  bcftools mpileup -d10000 -Ov -R $snps -o ${bam}.R.vcf -f $ref ${bam}.bam &
done

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
SL3.0ch00       244     .       C       A       10      PASS    AF=1
SL3.0ch00       247     .       C       A       10      PASS    AF=1

head -n 1000 $snps > snps_head.vcf
snps_head=snps_head.vcf
# vcftools --gzvcf $snps_head --freq --out chr1_analysis
# bcftools mpileup -d10000 -Ob -R $snps_head -o ${bam}.r.bcf -f $ref ${bam}.bam

# DP : combined depth across samples, e.g. DP=154
# PL : the phred-scaled genotype likelihoods rounded to the closest integer (and otherwise defined precisely asthe GL field) (Integers)
# In the vcf file, the I16 category is formatted as:
1 #reference Q13 bases on the forward strand
2 #reference Q13 bases on the reverse strand
3 #non-ref Q13 bases on the forward strand
4 #non-ref Q13 bases on the reverse strand
5 sum of reference base qualities
6 sum of squares of reference base qualities
7 sum of non-ref base qualities
8 sum of squares of non-ref base qualities
9 sum of ref mapping qualities
10 sum of squares of ref mapping qualities
11 sum of non-ref mapping qualities
12 sum of squares of non-ref mapping qualities
13 sum of tail distance for ref bases
14 sum of squares of tail distance for ref bases
15 sum of tail distance for non-ref bases
16 sum of squares of tail distance for non-ref 
!
install.packages("tidyverse")
install.packages("vcfR")
# HMM
install.packages("HMM")

library(vcfR)
library(HMM)
library(dplyr)
library(ggplot2)
library(tidyverse)
library("rtracklayer")

  # Let's see an example. There are two possible states called "Target" and "Outlier" in a test data, and their selecting probabilities are as below,

 mystates <- c("00","01", "11")
 mysymbols <- c("HZ_M2", "Hetero", "HZ_Penn")
 mytargetProb <- c(0.45, 0.1, 0.45)

 het2ho <- 1e-8
 het <- 1-2*het2ho
 ho2het <- 1e-9
 ho <- 1-2*ho2het
     # transProbs ‘transProbs’ is a (number of
     #      states)x(number of states)-sized matrix,
     #      which contains the transition probabilities
     #      between states.  The entry
     #      ‘transProbs[X,Y]’ gives the probability of
     #      a transition from state ‘X’ to state ‘Y’.
     #      The rows of the matrix must sum to 1.

 mytransProb <- matrix(c(ho, ho2het, ho2het,
                       het2ho, het, het2ho,
                       ho2het, ho2het, ho),
                     3,byrow=T)

 ehet <- 1/3
 eho <- 0.9
 eho2other <- (1-eho)/2

 myemissionP <- matrix(c(eho, eho2other , eho2other,
                       1/3, 1/3, 1/3,
                       eho2other, eho2other , eho),
                     3,byrow=T
                     )

 hmm <- HMM::initHMM(States = mystates, Symbols = mysymbols, startProbs=mytargetProb, transProbs=mytransProb, emissionProbs=myemissionP)
 tmp2$stateViterbi <- HMM::viterbi(hmm, tmp3$gt)


# Based on those selection probabilities, we build a transition probability matrix.

threshold_log <- 5 
basepath <-  "/projects/TRIBE/F4_genotyping/"
setwd(basepath)
# genomes <- c("Penn", "M82")

hetmaps <- list(mygenomes[1], mygenomes[2])
hetmaps_condensed <- list(SolLyc = list(), SolPen = list())
hetmaps_condensed2 <- list(SolLyc = list(), SolPen = list())
for (mygenome in mygenomes) {
  # mygenome <- mygenomes[1]
  subpath <- file.path("vcfs", mygenome)
  myfiles <- list.files(file.path(basepath, subpath),"R.vcf.gz$", full.names = F)

  for (myfile in myfiles){
    libname <- str_replace(myfile, "(.*)_Aligned.sortedByCoord.out.R.vcf.gz", "\\1")
    vcf <- read.vcfR( file.path(basepath, subpath, myfile), verbose = FALSE )

    I16 <- sapply(strsplit((vcf@fix[,8]),";"),function(x) x[2])
    I16b <- sapply(strsplit((I16),"="),function(x) x[2])
    ref <- sapply(strsplit((I16b),","),function(x) sum(as.numeric(x[1:2])))
    alt <- sapply(strsplit((I16b),","),function(x) sum(as.numeric(x[3:4])))

    PL1 <- as.integer(sapply(strsplit((vcf@gt[,2]),","),function(x) x[1]))
    PL2 <- as.integer(sapply(strsplit((vcf@gt[,2]),","),function(x) x[2]))
    PL3 <- as.integer(sapply(strsplit((vcf@gt[,2]),","),function(x) x[3]))

    tmp <- cbind(as.data.frame(vcf@fix), ref, alt, PL1, PL2, PL3)
    rm(ref, alt, PL1, PL2, PL3)

    tmp2 <- tmp %>%
      # tail() %>%
      # mutate(genotype2 =  pmap_int(list(PL1, PL2, PL3), function(...) which.max(c(...)))) %>%
      dplyr::filter(ref > threshold_log | alt > threshold_log) %>%
      dplyr::filter(!is.na(PL1)) %>%
      mutate(genotype =  pmap_int(list(PL1, PL2, PL3), function(...) which.min(c(...)))) %>%
      # mutate(mygt =  pmap_dbl(list(PL1, PL2, PL3), gt)) %>%
      mutate(log_ratio = log2((ref + 1) / (alt + 1))) %>%
      mutate(id = row_number()) %>%
      mutate(pos = as.integer(POS)) %>%
      mutate(gt = fct_recode(factor(genotype), "HZ_M82" = "1", "Hetero" = "2", "HZ_Penn" = "3"))

    # table(tmp2$genotype)

    # with(tmp2, plot(ref,alt, col = mygt, xlim = c(0,200), ylim = c(0,200)))
    tmp2$stateViterbi <- HMM::viterbi(hmm, tmp2$gt)
    hetmaps[[mygenome]][[libname]] <- tmp2
  }
}
saveRDS(hetmaps, file.path(basepath, file = "hetmaps.Rds"))


# count <- function(x) table(x$stateViterbi)/nrow(x)*100
# round(sapply(hetmaps[["SolLyc"]], count),1)
# round(sapply(hetmaps[["SolPen"]], count),1)


for (mygenome in mygenomes) {
  for (myfile in myfiles){
      # libname <- libnames[2]
    tmp2 <- hetmaps[[mygenome]][[libname]]
  gg <- ggplot(data = tmp2, aes(x = pos, y = log_ratio, color = stateViterbi, group = 1)) +
    geom_line() +
    facet_grid(CHROM~.,scale="free")
  ggsave(file=file.path(paste0(libname,"line_whichmin_viterbi.pdf")),gg,height=10,width=30)

    tmp2 <- hetmaps[[mygenome]][[libname]]

  tmp3 <- tmp2 %>%
    filter(CHROM == "SL3.0ch02" & as.integer(POS) < 4.055e7 & as.integer(POS) > 4.0544e7)
  gg <- ggplot(data = tmp3, aes(x = pos, y = log_ratio, color = gt, group = 1)) + geom_line()
p <- ggplotly(gg)
  library(plotly)
SL3.0ch02:40,544,663-40,544,702

  # xlim(4e7, 5e7)
  # facet_grid(CHROM~.,scale="free")
  ggsave(file=file.path(paste0(libname,"line_whichmin_chr6_viterbi2.pdf")),gg,height=3,width=30)


  gg <- ggplot(data = tmp2, aes(x = pos, y = log_ratio, color = gt, group = 1)) +
    geom_line() +
    facet_grid(CHROM~.,scale="free")
  ggsave(file=file.path(paste0(libname,"_line_whichmin2.pdf")),gg,height=10,width=30)
  rm(tmp,tmp2,tmp3,vcf,I16)
}

#------------------------------------ bin overlap of SNPs/genotypes
# load("/data/public_data/tomato/additional_resources/R_code_and_objects/R-objects-bins-SL30-Penn_v2_annot.rdata")
# contins bins for both slyc and penn
# mybins <- bins[["200"]]$gr
mygenomes <- c("SolLyc", "SolPen")
mygenome <- mygenomes[1]
# mygenomes <- c("M82", "Penne")
#------------------------------------ filtering genotypes
# filtering genotypes based on WT control
# The idea is to exclude all regions that were not corectly identified as parental genotype for the parent library
# This is done for each bin!
hetmaps <- readRDS(file = file.path(basepath, "hetmaps.Rds"))

#------------------------------------ reducing same adjacant SNPs in one interval
  # for (mygenome in mygenomes) {
    mygenome <- mygenomes[1]
    libnames <- names(hetmaps[[mygenome]])
    libnames <- libnames[-grep("DDM", libnames)]
    # libnames <- libnames[-grep("P1515", libnames)]
    libnames <- libnames[-grep("P1511", libnames)]

    for (libname in libnames){
      # libname <- libnames[1]
      dat <- hetmaps[[mygenome]][[libname]]
      # construct consecutive intervals (merge snps if adjacent ones are similar)

      dat2 <- dat %>%
        # head %>%
        mutate(seqnames = CHROM, 
               start = POS,
               end = POS,
               strand = "*",
               score = stateViterbi) %>%
        dplyr::select(seqnames, start, end, strand, score)

      gr <- makeGRangesFromDataFrame(dat2, keep.extra.columns=TRUE)
      gr2 <- gaps(gr) # find intervals between SNPs
      grre <- resize(gr, 2, fix = "end") # widen SNPs by 1 bp each side (no overlapping gaps above)
      grre2 <- resize(grre, 3, fix = "start") # widen SNPs by 1 bp each side (no overlapping gaps above)
      lst  <-  GenomicRanges::reduce(split(grre2, gr$score))
      gr2$score <- "None" 
      gr2$score[countOverlaps(gr2, lst$"00") == 2] <- "00"
      gr2$score[countOverlaps(gr2, lst$"01") == 2] <- "01"
      gr2$score[countOverlaps(gr2, lst$"11") == 2] <- "11"
      # Merging adjacent bins with same property in GenomicRanges object
      result <- unlist(GenomicRanges::reduce(split(c(gr2, grre2), ~score)))
      result$score <- names(result)
      result$ID <- paste(names(result), seq_along(result), sep = "_")
      result <- sort(subset(result, score != "None"))
      table(gr2$score, useNA = "ifany")
      hetmaps_condensed[[mygenome]][[libname]] <- result
      export(sort(result), paste0(libname,".gff"), format = "GFF3")
    }

    for (libname in libnames){
      message(libname)
      result <- hetmaps_condensed[[mygenome]][[libname]]
       print(round(100*prop.table(sapply(width(split(result, ~score)),sum)),1))
      # 00        01        11
      # 539402424 108883806 310027624
    }
saveRDS(hetmaps_condensed, file.path(basepath, file = "hetmaps_condensed.Rds"))

# condensing further: combining replicates only keeping positions that agree
hetmaps_condensed <- readRDS(file.path(basepath, file = "hetmaps_condensed.Rds"))

# dput(F4s)
plants <-  c("M82C", "PenneC", "P1512", "P2561", "P2562", "P3611", "P3612", "P4041", "P4042")
for (plant in plants) {
  # plant=plants[1]
  idx <- grep(plant, names(hetmaps_condensed[[mygenome]]))
  tmp <- hetmaps_condensed[[mygenome]][idx]
  tmp1 <- tmp[[1]]
  tmp2 <- tmp[[2]]
  ist <- intersect(split(tmp1, tmp1$score),
                   split(tmp2, tmp2$score))
  ist_stack <- stack(ist, index.var = "score")
  hetmaps_condensed2[[mygenome]][[plant]] <- ist_stack
  ist_stack$ID <- paste(ist_stack$score, seq_along(ist_stack$score), sep = "_")
  export(sort(ist_stack), paste0(plant, "_", mygenome, "_genotype.gff"), format = "GFF3")
}

saveRDS(hetmaps_condensed2, file.path(basepath, file = "hetmaps_condensed2.rds"))
