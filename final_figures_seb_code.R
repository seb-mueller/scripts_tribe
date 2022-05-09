
path_base <-  "/projects/TRIBE"
# laptop
path_base <-  "/home/seb/workspace/tribe"
source(file.path(path_base, "scripts/project_settings.R"))
# load(file=file.path(path_base, "current.workspace.Rdata"))
load(file = file.path(path_comparison, paste0(mygenome, "_slydf_plus_statsv39.rdata")))
load(file = file.path(path_comparison, paste0("R-objects-bins-SL30-Penn_v7_annot.rdata")))
# save(list = ls(), file=file.path(path_base, "current.workspace.Rdata"))
factor3     <-  "F4"
    #     factor2="TEs_Order"
    # factor2s <- c("TEs_Order", "TEs_Superfamily", "srnas_2124ratio_class", "isPromotor", "isTE", "annotation", "chromatin_state", "alltrue")
# factor2s <- c("TEs_Order_renamed", "dcl2_class2", "TEs_Order", "TEs_Superfamily", "srnas_2124ratio_class", "srnas_2124ratio_class", "annotation", "chromatin_state", "isLTRpred", "alltrue")
factor2s <- c("TEs_Order_renamed", "srnas_dcl2_class2", "rnaseq_dcl2_class2_promoter", "TEs_Superfamily", "rnaseq_dcl2_class2", "srnas_dcl2_class2", "srnas_2124ratio_class", "annotation", "srnas_2124ratio_class", "chromatin_state", "alltrue")
factor2s <- c("TEs_Order_renamed"

path_wd <- file.path(path_comparison, "sRNA_genomewide")
setwd(path_wd)
knitr::opts_knit$set(root.dir = path_wd)
# setwd("/home/sm934/workspace/tribe/scripts_tribe/plots/")
# rmarkdown::render("/home/sm934/workspace/tribe/scripts_tribe/general_slyc_stats.rmd")
# rmarkdown::render("/projects/TRIBE/scripts/general_slyc_stats.Rmd", output_dir = path_wd, output_file = paste("DRM", "general.html", sep = "_"))

#------------------- TEs non bins {{{
factor1 <- "class2"
factor2 <- "Superfamily"
factor2 <- "Order"

# factor2s <- c("Superfamily", "Order")

setwd(file.path(path_base_rnaseq, "20200217_RNAseq_parents_and_F4_seb_analyis_sara_mapped"))

rna_cnts_list <- list()
prefixes <- c("RNAseq_TEs_", "RNAseq_genes_")
prefix <- prefixes[1] # manual set til now..
prefix <- prefixes[2] # manual set til now..
rna_cnts_list[[prefix]] <- readRDS(file.path(path_base_rnaseq, "20200217_RNAseq_parents_and_F4_seb_analyis_sara_mapped", 
                           paste0(prefix , "SolLyc_edgeR_Robject.rds")))



# list_numbers <- list()
# for (myF4 in F4s_rnaseq) {
# myF4 <- F4s_rnaseq[1]
#   lrt <- rna_cnts_list[[prefix]][[myF4]]
# }


# this is using the non bin approach
list_perc <- list()
# only works for TEs!!
for (myF4 in F4s_rnaseq) {
  #   myF4 <- F4s_rnaseq[1]
  lrt <- rna_cnts_list[[prefix]][[myF4]]
  for (factor2 in factor2s) {

    mycounts <- lrt$table_all %>%
      dplyr::count(!!sym(factor1), !!sym(factor2)) %>%
      group_by(!!sym(factor1)) %>%
      mutate(!!sym(paste0(factor1, "_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      group_by(!!sym(factor2)) %>%
      mutate(!!sym(paste0(factor2,"_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      mutate(F4 = myF4)

    list_perc[[factor2]][[myF4]] <- mycounts

    write_csv(table_all, path = paste0(prefix, "stats_", myF4, "_", factor1,factor2, ".csv"))
    write_csv(mycounts , path = paste0(prefix, "stats_", myF4, "_", factor1,factor2, "_freq.csv"))
    # mybarplot2(factor1,factor2, mycounts = mycounts, mywidth = 16)

    gg <- ggplot(mycounts, aes(x = !!sym(factor1), 
                               y = !!sym(paste0(factor1, "_freq")), 
                               fill = fct_rev(!!sym(factor1)))) +
              #           gglayers +
              facet_grid(vars(), vars(!!sym(factor2))) +
              geom_bar(stat="identity", position="dodge") + 
              scale_fill_viridis_d() +
              theme_bw() +
              theme(legend.position = "right",
                    axis.ticks = element_blank(),
                    axis.text.x = element_text(angle = 300, hjust = 0))
              mytitle <- paste0(prefix, "stats_", myF4, "_rev2_", factor1,factor2, ".pdf")
              ggsave(mytitle, gg, height = 6, width = 9)

              gg <- ggplot(mycounts, aes(x = !!sym(factor2), 
                                         y = !!sym(paste0(factor2, "_freq")), 
                                         fill = fct_rev(!!sym(factor1)))) +
              gglayers +
              #           facet_grid(vars(), vars(!!sym(factor1))) +
              geom_bar(stat="identity", position="fill") + 
              scale_fill_viridis_d() +
              theme_bw() +
              theme(legend.position = "right",
                    axis.ticks = element_blank(),
                    axis.text.x = element_text(angle = 300, hjust = 0))
              mytitle <- paste0(prefix, "stats_", myF4, "_rev_", factor1,factor2, ".pdf")
              ggsave(mytitle, gg, height = 6, width = 4)
  }
}

# compiling all TE fractions in a data frame (for all F4s)
factor2 <- "Superfamily"
list_perc[[factor2]][[myF4]] <- mycounts
TEs_perc_df <- reduce(list_perc[[factor2]], rbind)
#   class2 Superfamily                                 n class2_freq Superfamily_freq F4
#   <chr>  <chr>                                   <int>       <dbl>            <dbl> <chr>
# 1 -DEG   CACTA                                       6        2.21            0.7   P4042
# 6 -DEG   EPRV                                        5        1.85            1.12  P4042


# EPRV plot


# those figures were not using bins, but lrt object!

TEs_perc_df_filter <- TEs_perc_df %>%
  #   filter(Superfamily == "EPRV") %>%
  filter(class2 != "inbetween")

    gg <- ggplot(TEs_perc_df_filter, aes(x = !!sym(factor1), 
                           y = !!sym(paste0(factor1, "_freq")))) +
#           gglayers +
          facet_grid(vars(), vars(!!sym(factor2))) +
          geom_boxplot() +
          geom_point(aes(color = !!sym(factor3))) +
                 theme_bw() +
                 theme(legend.position = "right",
                       axis.ticks = element_blank(),
                       axis.text.x = element_text(angle = 300, hjust = 0))

    mytitle <- paste0("TEs_stats_", "_allF4s_", factor1,factor2, ".pdf")
    ggsave(mytitle, gg, height = 6, width = 9)

# trying to create Superfamily plot for DESL for all F4 
# similar to TEs stats rev2
# bits and pieces:
#--- TEs non bins }}}

    myF4 ="P4042"
    infix = ""
    # gtfilter = "11" # for penn homozygous
    # gtfilter = "00" # for lyc
    # myprefix = "srnas_"
    # mydf = slydf
    myprefix <- "srnas_unique_"
    myprefix <- "srnas_"
    suffix <- "_class2"
    infix = "_PenneC"
    factor2  <- paste0(myprefix, myF4, infix, suffix)

    factor1="srnas_2124ratio_class"




# collecting all srna counts for each F4 as basis for plotting later
#
#------------------- bin approach for all data !!!
list_perc_srna2 <- list()
for (myF4 in F4s_srnas) {
  #   myF4 <- F4s_srnas[7]
  message(myF4)
  genotype <- sym(paste0("genotype_", myF4))
  #   rnaseq_a  <- paste0(mydata, myF4, "_class")
  #save number of DEG/DET for each F4
  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter)
  message(nrow(mydf_filtera))

  for (factor2 in factor2s) {
    #     factor2="TEs_Order_renamed"
    mycounts <- mydf_filtera %>%
      dplyr::count(!!sym(factor2)) %>%
      group_by(!!sym(factor2)) %>%
      mutate(!!sym(paste0(factor2,"_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      mutate(F4 = myF4)
    # remove F4 name to be able to merge all F4s in one df
    colnames(mycounts) <- str_replace(colnames(mycounts), paste0("_", myF4), "")

    list_perc_srna2[[factor2]][[myF4]] <- mycounts
  }
}
# save(list_perc_srna2, gtfilter, file="penn_list_perc_srnas.Rdata")

# plotting numbers of homozygous bins for all subcategories
# factor2="TEs_Order_renamed"
for (factor2 in factor2s) {
    perc_df_all2 <- reduce(list_perc_srna2[[factor2]], rbind) %>% 
      drop_na() %>%
      filter(!!sym(factor2) != "none")
    gg <- ggplot(perc_df_all2,
                # filter our only DEG/DESL/DMR
                aes(x = !!sym(factor2), 
                     y = n)) +
                scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3, digits = 2),
                                   breaks = scales::pretty_breaks(n = 8)) +
                geom_boxplot(alpha = 0.7) +
                geom_point(aes(shape = !!sym(factor3))) +
                scale_shape_manual(values = F4_shapes) +
                scale_fill_manual(values = TRIBE_colors) +
                theme_bw() +
                theme(legend.position = "right",
                      axis.ticks = element_blank(),
                      axis.text.x = element_text(angle = 300, hjust = 0))
                mytitle <- paste0("stats_", "_allF4s_", factor2, ".pdf")
                ggsave(filename=mytitle,plot= gg, height = 3, width = 3)
                #                 write_csv(perc_df_all2, path = paste0(mytitle, ".csv"))
}
#                   miRNAs_Slyc
# srnas_P4042_class2   none sly-miR10529 sly-MIR10529 sly-miR10530
#          -DESL       3037            0            0            0
#          +DESL      13833            0            0            1
#          inbetween 401817            0            0            0
#          non-DESL  961118            1            1            0
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-MIR10530 sly-miR10535a sly-miR10535b sly-MIR10535b
#          -DESL                0             1             0             0
#          +DESL                0             0             0             0
#          inbetween            0             0             1             1
#          non-DESL             1             0             0             0
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR10537 sly-MIR10537 sly-miR10538 sly-miR10539
#          -DESL                0            0            0            0
#          +DESL                0            0            0            2
#          inbetween            0            0            0            0
#          non-DESL             1            1            1            0
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR10540 sly-MIR10540 sly-MIR156d sly-miR156d-5p
#          -DESL                0            0           0              0
#          +DESL                0            0           0              1
#          inbetween            1            1           1              0
#          non-DESL             0            0           0              0
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR156e-3p sly-miR156e-5p sly-miR159b
#          -DESL                  0              0           0
#          +DESL                  0              1           0
#          inbetween              0              0           1
#          non-DESL               1              0           1
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR164a-3p sly-miR164b-3p sly-miR166c-3p
#          -DESL                  0              0              0
#          +DESL                  0              0              0
#          inbetween              0              1              0
#          non-DESL               1              0              1
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR166c-5p sly-MIR167b sly-miR167b-3p
#          -DESL                  0           0              0
#          +DESL                  0           0              0
#          inbetween              0           0              1
#          non-DESL               1           1              0
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR167b-5p sly-miR168a-3p sly-miR168a-5p
#          -DESL                  0              0              0
#          +DESL                  0              0              0
#          inbetween              1              1              1
#          non-DESL               0              0              0
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-MIR169e sly-miR169e-3p sly-miR171e sly-miR171f
#          -DESL               0              1           0           0
#          +DESL               0              0           1           0
#          inbetween           0              0           0           0
#          non-DESL            1              0           0           1
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR172d sly-MIR172d sly-miR1916 sly-MIR1916
#          -DESL               0           0           0           0
#          +DESL               0           0           0           0
#          inbetween           0           0           0           0
#          non-DESL            1           1           1           1
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR1917 sly-miR319b sly-MIR319b sly-miR319d
#          -DESL               0           0           0           0
#          +DESL               0           0           0           0
#          inbetween           1           0           0           0
#          non-DESL            0           1           1           1
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-MIR319d sly-MIR394 sly-miR394-3p sly-miR408
#          -DESL               0          0             0          0
#          +DESL               0          0             0          1
#          inbetween           1          1             1          0
#          non-DESL            0          0             0          0
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-MIR408 sly-miR482a sly-miR530 sly-MIR530
#          -DESL              0           0          0          0
#          +DESL              0           0          0          0
#          inbetween          1           0          1          1
#          non-DESL           0           1          0          0
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR5304 sly-MIR5304 sly-miR6022 sly-MIR6022
#          -DESL               0           0           0           0
#          +DESL               0           0           1           0
#          inbetween           0           1           0           1
#          non-DESL            1           0           0           0
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR7981c sly-MIR7981c sly-miR7981e sly-MIR7981e
#          -DESL                0            0            0            0
#          +DESL                0            0            0            0
#          inbetween            1            0            0            1
#          non-DESL             0            1            1            0
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR7981f sly-MIR7981f sly-miR9469-3p
#          -DESL                0            0              0
#          +DESL                0            0              0
#          inbetween            0            1              0
#          non-DESL             1            0              1
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR9469-5p sly-MIR9470 sly-miR9470-3p
#          -DESL                  0           0              0
#          +DESL                  0           0              0
#          inbetween              1           1              0
#          non-DESL               0           1              1
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR9470-5p sly-miR9472-3p sly-miR9472-5p
#          -DESL                  0              0              0
#          +DESL                  0              1              1
#          inbetween              0              0              0
#          non-DESL               1              0              0
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR9474-3p sly-miR9474-5p sly-miR9476-5p
#          -DESL                  0              0              0
#          +DESL                  0              0              0
#          inbetween              1              0              1
#          non-DESL               0              1              0
#                   miRNAs_Slyc
# srnas_P4042_class2 sly-miR9479-3p sly-miR9479-5p
#          -DESL                  0              0
#          +DESL                  0              0
#          inbetween              0              0
#          non-DESL               1              1

slydf %>%
  # filter(miRNAs_Slyc == "sly-miR482a")
  filter(miRNAs_Slyc == "sly-miR6026")

with(mydf_filtera, table(srnas_P4042_class2, miRNAs_Slyc))
with(mydf_filtera, table(miRNAs_Slyc, srnas_P4042_class2))
with(mydf_filtera, table(srnas_P4042_class2, isRN))

with(mydf_filtera, table(srnas_unique_P4042_class2))
# srnas_unique_P4042_class2
#     -DESL     +DESL inbetween  non-DESL 
#        23      1013      1041     10310 
with(slydf, table(srnas_unique_P4042_class2, genotype_P4042))
with(slydf, table(srnas_P4042_class2, genotype_P4042))
with(slydf, table(srnas_P4042_class2, genotype_P4042, TEs_Order_renamed))
# , , TEs_Order_renamed = Pararetrovirus
#                   genotype_P4042
# srnas_P4042_class2     00     01     11   none
#          -DESL        189     64    505    112
#          +DESL       1408    340   1777    354
#          inbetween   4068   1011   5628    828
#          non-DESL    4212   1152   6063    742
with(slydf, table(srnas_unique_P4042_class2, genotype_P4042, TEs_Order_renamed))
                         # genotype_P4042
# srnas_unique_P4042_class2     00     01     11   none
                # -DESL         30     18    335     43
                # +DESL        125     17     50     11
                # inbetween    289     65    434     67
                # non-DESL    9433   2467  13154   1915

with(mydf_filtera, table(srnas_unique_P4042_class2, TEs_Order_renamed))
with(mydf_filtera, table(srnas_P4042_class2, TEs_Order_renamed, chromatin_state))
#                          TEs_Order_renamed
# srnas_unique_P4042_class2  TIR  LTR none LINE Pararetrovirus Helitron Other_TEs
#                 -DESL        0    2   11    0              4        0         4
#                 +DESL       52  144  554   19             47       20       148
#                 inbetween   64  191  485   30             38        7       166
#                 non-DESL   528 3662 3142  119            552       43      1344
with(mydf_filtera, table(srnas_P4042_class2))
# srnas_P4042_class2
#     -DESL     +DESL inbetween  non-DESL 
#        59      4350      5149      2829 
with(mydf_filtera, table(srnas_P4042_class2, TEs_Order_renamed))
#                   TEs_Order_renamed
# srnas_P4042_class2  TIR  LTR none LINE Pararetrovirus Helitron Other_TEs
#          -DESL        0   10   23    1             11        0        12
#          +DESL      229 1676 1003   33            403       42       487
#          inbetween  287 1473 2002   98            166       20       765
#          non-DESL   128  840 1164   36             61        8       398
with(mydf_filtera, table(srnas_P4042_class2, TEs_rnaseq_P4042_class2))
#                   TEs_rnaseq_P4042_class2
# srnas_P4042_class2 -DET +DET inbetween non-DET none
#          -DESL        0    0         0       7   52
#          +DESL        6   23       228     847 3246
#          inbetween    0    3       126     733 4287
#          non-DESL     1    0        60     412 2356

mydf_filtera %>%
  filter(srnas_P4042_class2 == "+DESL") %>%
  filter(TEs_rnaseq_P4042_class2 == "+DET") %>%
  select(c(1:3,11:14))

with(mydf_filtera, table(TEs_rnaseq_P4042_class2))
# TEs_rnaseq_P4042_class2
#      -DET      +DET inbetween   non-DET      none 
#         7        26       414      1999      9941 
with(mydf_filtera, table(srnas_P3612_class2, TEs_rnaseq_P3612_class2))

# collecting all srna counts for each F4 as basis for plotting later
#------------------- collecting all srna counts for each F4 as basis for plotting later {{{
    myF4 ="P4042"
    infix = ""
    # gtfilter = "11" # for penn homozygous
    # gtfilter = "00" # for lyc
    # myprefix = "srnas_"
    # mydf = slydf
    myprefix <- "srnas_unique_"
    myprefix <- "srnas_"
    suffix <- "_class2"
    infix = "_PenneC"
    factor2  <- paste0(myprefix, myF4, infix, suffix)

list_perc_srna <- list()
list_cnt_srna <- list()
for (myF4 in F4s_srnas) {
  #   myF4 <- F4s_rnaseq[1]
  message(myF4)
  genotype <- sym(paste0("genotype_", myF4))
  factor1  <- paste0(myprefix, myF4, infix, suffix)
  #   rnaseq_a  <- paste0(mydata, myF4, "_class")
  #save number of DEG/DET for each F4
  # mydf_filtera <- slydf %>%
    # filter(!!genotype == gtfilter)
  # message(nrow(mydf_filtera))

  mycounts2 <- slydf %>%
    filter(!!genotype == gtfilter) %>%
      dplyr::count(!!sym(factor1)) %>%
      mutate(F4 = myF4)
  colnames(mycounts2) <- str_replace(colnames(mycounts2), paste0("_", myF4), "")
 list_cnt_srna[[myF4]] <- mycounts2
}

  for (factor2 in factor2s) {
    #     factor2="TEs_Order_renamed"
    #     factor2="rnaseq_dcl2_class2_promoter"
    mycounts <- mydf_filtera %>%
      dplyr::count(!!sym(factor1), !!sym(factor2)) %>%
      drop_na() %>%
      filter(!!sym(factor2) != "none") %>%
      group_by(!!sym(factor1)) %>%
      mutate(!!sym(paste0(factor1, "_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      group_by(!!sym(factor2)) %>%
      mutate(!!sym(paste0(factor2,"_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      mutate(F4 = myF4)
    # remove F4 name to be able to merge all F4s in one df
    colnames(mycounts) <- str_replace(colnames(mycounts), paste0("_", myF4), "")

    list_perc_srna[[factor2]][[myF4]] <- mycounts
  }
}
save(list_perc_srna2, list_perc_srna, gtfilter, file="penn_list_perc_srnas.Rdata")
load(file.path(path_base, "comparisons_using_big_table/sRNA_genomewide/penn/penn_list_perc_srnas.Rdata"))
# MANUAL: do the same for list_perc_srna_uniqe
# setting myprefix="srnas_unique_"
# list_perc_srna_unique <- list_perc_srna
# list_cnt_srna_unique <- list_cnt_srna
#--- collecting all srna counts for each F4 as basis for plotting later }}}

#------------------- DMRs bins {{{
setwd(file.path(path_comparison, "DMR_genomewide"))
# collecting all DMR counts for each F4 as basis for plotting later
    myprefix <- "meth_"
    suffix <- "_class_stringent"
    infix <- "CHH"

list_perc_DMR_CHH <- list()
list_cnt_DMR_CHH <- list()
list_perc_DMR_CHG <- list()
list_cnt_DMR_CHG <- list()
list_perc_DMR_CpG <- list()
list_cnt_DMR_CpG <- list()

for (myF4 in F4s_meth2) {
  #   myF4 <- F4s_meth[1]
  message(myF4)
  genotype <- sym(paste0("genotype_", myF4))
  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter)
  message(nrow(mydf_filtera))

  for (infix in contexts) {
    # infix = "CpG"
    factor1  <- paste0(myprefix, myF4, "_", infix, suffix)
    #   rnaseq_a  <- paste0(mydata, myF4, "_class")
    #save number of DEG/DET for each F4

    mycounts2 <- mydf_filtera %>%
      dplyr::count(!!sym(factor1)) %>%
      mutate(F4 = myF4)
    colnames(mycounts2) <- str_replace(colnames(mycounts2), paste0("_", myF4), "")
    #     list_cnt_DMR[[myF4]] <- mycounts2
      switch(infix,
             CHH = list_cnt_DMR_CHH[[myF4]] <- mycounts2,
             CHG = list_cnt_DMR_CHG[[myF4]] <- mycounts2,
             CpG = list_cnt_DMR_CpG[[myF4]] <- mycounts2)

    for (factor2 in factor2s) {
      #     factor2="TEs_Order_renamed"

      mycounts <- mydf_filtera %>%
        dplyr::count(!!sym(factor1), !!sym(factor2)) %>%
        drop_na() %>%
        filter(!!sym(factor2) != "none") %>%
        group_by(!!sym(factor1)) %>%
        mutate(!!sym(paste0(factor1, "_freq")) := round(100 * n / sum(n),2)) %>%
        ungroup() %>%
        group_by(!!sym(factor2)) %>%
        mutate(!!sym(paste0(factor2,"_freq")) := round(100 * n / sum(n),2)) %>%
        ungroup() %>%
        mutate(F4 = myF4)
      # remove F4 name to be able to merge all F4s in one df
      colnames(mycounts) <- str_replace(colnames(mycounts), paste0("_", myF4), "")

      switch(infix,
             CHH = list_perc_DMR_CHH[[factor2]][[myF4]] <- mycounts,
             CHG = list_perc_DMR_CHG[[factor2]][[myF4]] <- mycounts,
             CpG = list_perc_DMR_CpG[[factor2]][[myF4]] <- mycounts)

    }
  }
}

#--- DMRs bins }}}

#------------------- same as above (DESL) but for DEG and DETs {{{

    # mydf = slydf
    myprefixdet <- "TEs_rnaseq_"
    myprefixdeg <- "rnaseq_"
    suffix <- "_class2"
    myF4 ="P4042"
    # gtfilter = "11" # for penn homozygous
    # gtfilter = "00" # for lyc
    infix = ""
    infix = "_PenneC"

list_perc_dets <- list()
list_cnt_dets <- list()
list_perc_degs <- list()
list_cnt_degs <- list()
for (myF4 in F4s_rnaseq) {
  #   myF4 <- F4s_rnaseq[1]

  message(myF4)
  genotype <- sym(paste0("genotype_", myF4))
  factor1det  <- paste0(myprefixdet, myF4, infix, suffix)
  factor1deg  <- paste0(myprefixdeg, myF4, infix, suffix)
# * Column `TEs_rnaseq_P1512_PenneC_class2` is not found.
  #   rnaseq_a  <- paste0(mydata, myF4, "_class")
  #save number of DEG/DET for each F4
  # mydf_filtera <- slydf %>%
  #   filter(!!genotype == gtfilter)
  # message(nrow(mydf_filtera))

  # for troubleshooting
  # mycounts2distinct <- mydf_filtera %>%
  # mycounts2<- mydf_filtera %>%
  #     distinct(TEs_name, .keep_all= TRUE) %>%
  #     dplyr::count(!!sym(factor1det)) %>%
  #     mutate(F4 = myF4)
  # TEs_rnaseq_P4042_class2      n    F4
# 1                    -DET     54 P4042
# 2                    +DET    208 P4042
# 3               inbetween   1649 P4042
# 4                 non-DET  54990 P4042
# 5                    none 139687 P4042

  colnames(mycounts2) <- str_replace(colnames(mycounts2), paste0("_", myF4), "")

 list_cnt_dets[[myF4]] <- mycounts2
  mycounts2 <- slydf %>%
      filter(!!genotype == gtfilter) %>%
      dplyr::count(!!sym(factor1det)) %>%
      mutate(F4 = myF4)
  colnames(mycounts2) <- str_replace(colnames(mycounts2), paste0("_", myF4), "")
 list_cnt_dets[[myF4]] <- mycounts2

  mycounts2 <- slydf %>%
      filter(!!genotype == gtfilter) %>%
      distinct(Genes_Slyc, .keep_all= TRUE) %>%
      dplyr::count(!!sym(factor1deg)) %>%
      mutate(F4 = myF4)
  colnames(mycounts2) <- str_replace(colnames(mycounts2), paste0("_", myF4), "")

 list_cnt_degs[[myF4]] <- mycounts2
}

  # for DETs
  for (factor2 in factor2s) {
    #     factor2="TEs_Order_renamed"
    mycounts <- mydf_filtera %>%
      dplyr::count(!!sym(factor1det), !!sym(factor2)) %>%
      group_by(!!sym(factor1det)) %>%
      mutate(!!sym(paste0(factor1det, "_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      group_by(!!sym(factor2)) %>%
      mutate(!!sym(paste0(factor2,"_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      mutate(F4 = myF4)
    # remove F4 name to be able to merge all F4s in one df
    colnames(mycounts) <- str_replace(colnames(mycounts), paste0("_", myF4), "")

    list_perc_dets[[factor2]][[myF4]] <- mycounts
  }

  # for DEGs
  for (factor2 in factor2s) {
    #     factor2="TEs_Order"
    mycounts <- mydf_filtera %>%
      dplyr::count(!!sym(factor1deg), !!sym(factor2)) %>%
      group_by(!!sym(factor1deg)) %>%
      mutate(!!sym(paste0(factor1deg, "_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      group_by(!!sym(factor2)) %>%
      mutate(!!sym(paste0(factor2,"_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      mutate(F4 = myF4)
    # remove F4 name to be able to merge all F4s in one df
    colnames(mycounts) <- str_replace(colnames(mycounts), paste0("_", myF4), "")

    list_perc_degs[[factor2]][[myF4]] <- mycounts
  }
}
#--- same as above (DESL) but for DEG and DETs }}}
# save(list_perc_srna2, list_perc_srna, gtfilter, file="penn_list_perc_srnas.Rdata")

#------------------------------------2D: counts for pairwise combination for all F4s {{{

myprefix <- "srnas_unique_" # do this manually
myprefix <- "srnas_"
infix = ""
suffix <- "_class2"

list_perc_srna_2d <- list()
list_cnt_srna_2d <- list()
for (myF4 in F4s_srnas) {
  #   myF4 <- F4s_srnas[1]
  message(myF4)
  genotype <- sym(paste0("genotype_", myF4))
  factor1  <- paste0(myprefix, myF4, infix, suffix)
  #   rnaseq_a  <- paste0(mydata, myF4, "_class")
  #save number of DEG/DET for each F4
  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter)
  message(nrow(mydf_filtera))

  mycounts2 <- mydf_filtera %>%
      dplyr::count(!!sym(factor1)) %>%
      mutate(F4 = myF4)
  colnames(mycounts2) <- str_replace(colnames(mycounts2), paste0("_", myF4), "")
 list_cnt_srna_2d[[myF4]] <- mycounts2

  for (factor2a in factor2s) {
    for (factor2b in factor2s) {
      #     factor2a=factor2s[1]
      #     factor2b=factor2s[3]
      mycounts <- mydf_filtera %>%
        dplyr::count(!!sym(factor1), !!sym(factor2a), !!sym(factor2b)) %>%
        group_by(!!sym(factor1)) %>%
        mutate(!!sym(paste0(factor1, "_freq")) := round(100 * n / sum(n),2)) %>%
        ungroup() %>%
        group_by(!!sym(factor2a)) %>%
        mutate(!!sym(paste0(factor2a,"_freq")) := round(100 * n / sum(n),2)) %>%
        ungroup() %>%
        group_by(!!sym(factor2b)) %>%
        mutate(!!sym(paste0(factor2b,"_freq")) := round(100 * n / sum(n),2)) %>%
        ungroup() %>%
        mutate(F4 = myF4)
      # remove F4 name to be able to merge all F4s in one df
      colnames(mycounts) <- str_replace(colnames(mycounts), paste0("_", myF4), "")

      list_perc_srna_2d[[factor2a]][[factor2b]][[myF4]] <- mycounts
    }
  }
}

list_perc_srna_2d2 <- list()
list_cnt_srna_2d2 <- list()
for (myF4 in F4s_srnas) {
  #   myF4 <- F4s_srnas[1]
  message(myF4)
  genotype <- sym(paste0("genotype_", myF4))
  factor1  <- paste0(myprefix, myF4, infix, suffix)
  #   rnaseq_a  <- paste0(mydata, myF4, "_class")
  #save number of DEG/DET for each F4
  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter)
  message(nrow(mydf_filtera))

  mycounts2 <- mydf_filtera %>%
      dplyr::count(!!sym(factor1)) %>%
      mutate(F4 = myF4)
  colnames(mycounts2) <- str_replace(colnames(mycounts2), paste0("_", myF4), "")
 list_cnt_srna_2d[[myF4]] <- mycounts2

  for (factor2a in factor2s) {
    for (factor2b in factor2s) {
      #     factor2a=factor2s[1]
      #     factor2b=factor2s[3]
      mycounts <- mydf_filtera %>%
        dplyr::count(!!sym(factor1), !!sym(factor2a), !!sym(factor2b)) %>%
        group_by(!!sym(factor1)) %>%
        mutate(!!sym(paste0(factor1, "_freq")) := round(100 * n / sum(n),2)) %>%
        ungroup() %>%
        group_by(!!sym(factor2a)) %>%
        mutate(!!sym(paste0(factor2a,"_freq")) := round(100 * n / sum(n),2)) %>%
        ungroup() %>%
        group_by(!!sym(factor2b)) %>%
        mutate(!!sym(paste0(factor2b,"_freq")) := round(100 * n / sum(n),2)) %>%
        ungroup() %>%
        mutate(F4 = myF4)
      # remove F4 name to be able to merge all F4s in one df
      colnames(mycounts) <- str_replace(colnames(mycounts), paste0("_", myF4), "")

      list_perc_srna_2d[[factor2a]][[factor2b]][[myF4]] <- mycounts
    }
  }
}
# list_perc_srna_2d_unique  <- list_perc_srna_2d

# [1] "/projects/TRIBE/comparisons_using_big_table/sRNA_genomewide/another_run_nov20"
save(list=ls(pattern="*list*"), file="final_figures_seb_code2d.Rdata")
#}}}
# save(list=ls(pattern="*list*"), file="penn_final_figures_seb_code2d.Rdata")



# /* list_deg_table <- map(list_perc_srna[[factor2]], ~table(.x[,2])) */
# TEs_perc_df <- reduce(list_perc_srna[[factor2]],  bind(.x,.y))
# colnames(list_perc_srna[[factor2]][[1]])
# perc_df_all <- rbindlist(list_perc_srna[[factor2]])

#------------------------------------ plotting all F4s

#------------------------------------ first just number of DESLs

# this will produce barplots showing the number of DESL for each category and F4
load(file="final_figures_seb_code2d.Rdata")

factor1s <- c("srnas_class2", "srnas_unique_class2", "TEs_rnaseq_class2", "rnaseq_class2", "meth_CHH_class_stringent", "meth_CpG_class_stringent", "meth_CHG_class_stringent") 
factor1s <- c("srnas_class2", "TEs_rnaseq_class2", "rnaseq_class2") 
# factor1s <- c("srnas_class2", "TEs_rnaseq_class2", "rnaseq_class2") 
# factor1s <- c("srnas_unique_class2")
# factor1s <- c("srnas_class2")
# factor1s <- c("srnas_PenneC_class2")
cnt_df                                    <- list()
cnt_df2                                   <- list()
avg_de_list                               <- list()
avg_de_filter                             <- list()
perc_df_all_lists                         <- list()
perc_df_all_lists_2d                      <- list()
cnt_df[["srnas_class2"]]                  <- reduce(list_cnt_srna, rbind)
cnt_df[["srnas_PenneC_class2"]]                  <- reduce(list_cnt_srna, rbind)
cnt_df[["srnas_unique_class2"]]           <- reduce(list_cnt_srna_unique, rbind)
# cnt_df[["srnas_unique_class2"]]           <- reduce(list_cnt_srna, rbind)
cnt_df[["TEs_rnaseq_class2"]]             <- reduce(list_cnt_dets, rbind)
cnt_df[["rnaseq_class2"]]                 <- reduce(list_cnt_degs, rbind)
cnt_df[["TEs_rnaseq_PenneC_class2"]]             <- reduce(list_cnt_dets, rbind)
cnt_df[["rnaseq_PenneC_class2"]]                 <- reduce(list_cnt_degs, rbind)
cnt_df[["meth_CHH_class_stringent"]]      <- reduce(list_cnt_DMR_CHH, rbind)
cnt_df[["meth_CpG_class_stringent"]]      <- reduce(list_cnt_DMR_CpG, rbind)
cnt_df[["meth_CHG_class_stringent"]]      <- reduce(list_cnt_DMR_CHG, rbind)
avg_de_list[["srnas_class2"]]             <- reduce(list_perc_srna[["alltrue"]], rbind)
avg_de_list[["srnas_PenneC_class2"]]             <- reduce(list_perc_srna[["alltrue"]], rbind)
avg_de_list[["srnas_unique_class2"]]      <- reduce(list_perc_srna_unique[["alltrue"]], rbind)
avg_de_list[["TEs_rnaseq_class2"]]        <- reduce(list_perc_dets[["alltrue"]], rbind)
avg_de_list[["rnaseq_class2"]]            <- reduce(list_perc_degs[["alltrue"]], rbind)
avg_de_list[["TEs_PenneC_rnaseq_class2"]]        <- reduce(list_perc_dets[["alltrue"]], rbind)
avg_de_list[["rnaseq_PenneC_class2"]]            <- reduce(list_perc_degs[["alltrue"]], rbind)
avg_de_list[["meth_CHH_class_stringent"]] <- reduce(list_perc_DMR_CHH[["alltrue"]], rbind)
avg_de_list[["meth_CpG_class_stringent"]] <- reduce(list_perc_DMR_CpG[["alltrue"]], rbind)
avg_de_list[["meth_CHG_class_stringent"]] <- reduce(list_perc_DMR_CHG[["alltrue"]], rbind)

perc_df_all_lists[["srnas_class2"]]      <- list_perc_srna
perc_df_all_lists[["srnas_PenneC_class2"]]      <- list_perc_srna
perc_df_all_lists[["srnas_unique_class2"]]      <- list_perc_srna_unique
# perc_df_all_lists[["srnas_unique_class2"]]      <- list_perc_srna
perc_df_all_lists[["TEs_rnaseq_class2"]] <- list_perc_dets
perc_df_all_lists[["rnaseq_class2"]]     <- list_perc_degs
perc_df_all_lists[["meth_CHH_class_stringent"]] <- list_perc_DMR_CHH
perc_df_all_lists[["meth_CHG_class_stringent"]] <- list_perc_DMR_CHG
perc_df_all_lists[["meth_CpG_class_stringent"]] <- list_perc_DMR_CpG
perc_df_all_lists_2d[["srnas_class2"]] <- list_perc_srna_2d
perc_df_all_lists_2d[["srnas_PenneC_class2"]] <- list_perc_srna_2d
perc_df_all_lists_2d[["srnas_unique_class2"]] <- list_perc_srna_2d_unique
# perc_df_all_lists_2d[["TEs_rnaseq_class2"]] <- list_perc_dets
# perc_df_all_lists_2d[["rnaseq_class2"]] <- list_perc_degs
#2D
cnt_df_2d <- list()
tmp <- list_perc_srna_2d[[1]]
cnt_df_2d[["srnas_class2"]] <- reduce(list_perc_srna_2d[["alltrue"]][[1]], rbind)
cnt_df_2d[["srnas_class2"]] <- reduce(tmp, rbind)
cnt_df_2d[["srnas_unique_class2"]] <- reduce(list_perc_srna_2d[["alltrue"]][[1]], rbind)
cnt_df_2d[["srnas_unique_class2"]] <- reduce(tmp, rbind)


# working out the average fractions of DE elements as a baseline
# tmp: factor1s <- c("meth_CHH_class_stringent", "meth_CpG_class_stringent", "meth_CHG_class_stringent") 
# tmp: factor1s <- c("TEs_rnaseq_class2")
# tmp: factor1 <- c("srnas_class2")
# tmp: factor1 <- c("srnas_class2")
for(factor1 in factor1s) {
  # factor1 = factor1s[2]
  tmp <- avg_de_list[[factor1]] %>%
    group_by(!!sym(paste0(factor1))) %>%
    summarise(n = n(),
              mean = mean(!!sym(paste0("alltrue", "_freq")), na.rm = TRUE),
              median = median(!!sym(paste0("alltrue", "_freq")), na.rm = TRUE)) %>%
    filter(str_detect(!!sym(factor1), "^[+-].*")) %>%
    as.data.frame()
  tmp["ratio",3:4] <- tmp[1,3:4]/tmp[2,3:4]
  avg_de_filter[[factor1]] <- tmp
}

save(avg_de_filter, cnt_df, avg_de_list, perc_df_all_lists, file="final_figures_seb_code2dv2.Rdata")
load(file="final_figures_seb_code2dv2.Rdata")
load(file.path(path_base, "comparisons_using_big_table/sRNA_genomewide/penn/final_figures_seb_code2dv2.Rdata"))
load(file.path(path_base, "comparisons_using_big_table/sRNA_genomewide/penn/final_figures_seb_code2d.Rdata"))
load(file.path(path_base, "comparisons_using_big_table/sRNA_genomewide/penn/penn_final_figures_seb_code2d.Rdata"))
# $srnas_PenneC_class2
#       srnas_PenneC_class2  n      mean    median
# 1                   -DESL  7 0.2414286 0.2100000
# 2                   +DESL  7 0.4014286 0.3700000
# ratio                <NA> NA 0.6014235 0.5675676

#       meth_CHH_class_stringent n   mean median
# 1                         -DMR 6 0.0350   0.03
# 2                         +DMR 6 0.0800   0.06
# ratio                     <NA> 1 0.4375   0.50

# $srnas_class2
# # A tibble: 2 x 4
#   srnas_class2     n  mean median
#   <fct>        <int> <dbl>  <dbl>
# 1 -DESL            7 0.313   0.31
# 2 +DESL            7 1.16    1.16
# 3 NA               1 0.269  0.267
# 
# $TEs_rnaseq_class2
# # A tibble: 2 x 4
#   TEs_rnaseq_class2     n   mean median
#   <fct>             <int>  <dbl>  <dbl>
# 1 -DET                  7 0.0414   0.04
# 2 +DET                  7 0.11     0.11
# 
# $rnaseq_class2
# # A tibble: 2 x 4
#   rnaseq_class2     n  mean median
#   <fct>         <int> <dbl>  <dbl>
# 1 -DEG              7 0.541   0.51
# 2 +DEG              7 0.547   0.44
#   srnas_class2     n   mean median
#   <fct>        <int>  <dbl>  <dbl>
# 1 -DESL            7  0.313   0.31
# 2 +DESL            7  1.16    1.16
# 3 inbetween        7 24.1    29.1
# 4 non-DESL         7 74.4    69.6


# TODO:DEG should be made unique!?
for(factor1 in factor1s) {
  # factor1 = "rnaseq_class2"
  # factor1 = "rnaseq_P4042_PenneC_class2"
  # factor1 = "rnaseq_PenneC_class2"
  # factor1 = "TEs_rnaseq_PenneC_class2"
cnt_df2[[factor1]] <- cnt_df[[factor1]] %>%
  filter(!!sym(paste0(factor1)) != "inbetween") %>%
  filter(!!sym(paste0(factor1)) != "non-DESL") %>%
  filter(!!sym(paste0(factor1)) != "non-DET") %>%
  filter(!!sym(paste0(factor1)) != "non-DEG") %>%
  filter(!!sym(paste0(factor1)) != "non-DMR") %>%
  filter(!!sym(paste0(factor1)) != "none")
}


# plotting stacked barplots for fig_desl fig_deg, (first panels)
# for(factor1 in factor1s) {
for(factor1 in names(cnt_df)) {
  # factor1 <- factor1s[1]
  # factor1 <- names(cnt_df)[3]
# factor1tmp <- "TEs_rnaseq_class2"
# factor1tmp <- "_rnaseq_class2"
  mytitle <- paste0(factor1, "_number_using_bins_filtered.pdf")
  gg <- ggplot(
               cnt_df[[factor1]] %>%
               # filter(!!sym(paste0(factor1))!= "none")
               filter(TEs_rnaseq_class2 != "none")
             ,
               aes(x = F4, y = n, 
                   fill = fct_rev(!!sym(paste0(factor1))))) +
        gglayers2 +
        scale_fill_manual(values = TRIBE_colors) +
        scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3, digits = 2), breaks = scales::pretty_breaks(n = 8)) +
        geom_bar(stat="identity") +
        labs(title = mytitle)
  ggsave(mytitle, gg, height = 6, width = 5)
  write_csv(cnt_df[[factor1]], path = paste0(mytitle, ".csv"))

  # review: fig 1d/2a % instead of total numbers
  mytitle <- paste0(factor1, "_number_using_bins_filtered_pct.pdf")
  gg <- ggplot(
               tmp <- cnt_df[[factor1]] %>%
                 # filter(!!sym(paste0(factor1))!= "none") %>%
               filter(!!sym(factor1tmp) != "none") %>%
                 group_by(F4) %>%
                 mutate(perc=n/sum(n))
               ,
               aes(x = F4, y = perc, 
                   fill = fct_rev(!!sym(paste0(factor1))))) +
                   # fill = fct_rev(!!sym(factor1tmp)))) +
        gglayers2 +
        scale_fill_manual(values = TRIBE_colors) +
        scale_y_continuous(labels = scales::percent) +
        #scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3, digits = 2), breaks = scales::pretty_breaks(n = 8)) +
        geom_bar(stat="identity") +
        labs(title = mytitle)
  ggsave(mytitle, gg, height = 6, width = 4)
  write_csv(tmp, path = paste0(mytitle, ".csv"))

  mytitle <- paste0(factor1, "_number_using_bins_filtered_no-non.pdf")
  gg <- ggplot(cnt_df2[[factor1]], 
               aes(x = F4, y = n, 
                   fill = fct_rev(!!sym(paste0(factor1))))) +
                gglayers2 +
                scale_fill_manual(values = TRIBE_colors) +
                geom_bar(stat="identity", position = "dodge") +
                labs(title = mytitle)
  ggsave(mytitle, gg, height = 6, width = 5) #w =4 for non DMR

  mytitle <- paste0(factor1, "_number_using_bins_filtered_no-non_pct.pdf")
  gg <- ggplot(
               tmp <- cnt_df[[factor1]] %>%
                 group_by(F4) %>%
                 mutate(perc=n/sum(n)) %>%
                 filter(!!sym(paste0(factor1)) != "inbetween") %>%
                 filter(!!sym(paste0(factor1)) != "non-DESL") %>%
                 filter(!!sym(paste0(factor1)) != "non-DET") %>%
                 filter(!!sym(paste0(factor1)) != "non-DEG") %>%
                 filter(!!sym(paste0(factor1)) != "non-DMR") %>%
                 filter(!!sym(paste0(factor1)) != "none")
               ,
               aes(x = F4, y = perc, 
                   fill = fct_rev(!!sym(paste0(factor1))))) +
                   # fill = fct_rev(!!sym(factor1tmp)))) +

              gglayers2 +
              scale_fill_manual(values = TRIBE_colors) +
              scale_y_continuous(labels = scales::percent) +
              geom_bar(stat="identity", position = "dodge") +
              labs(title = mytitle)
            ggsave(mytitle, gg, height = 6, width = 4) #w =4 for non DMR
  write_csv(tmp, path = paste0(mytitle, ".csv"))

}

factor2 <- "alltrue"

    # those are plotting percentages per subclass as sused in fig_DEG and fig_DESL
    # i.e. very important. Sometimes scale and width are not optimal and need to be set
    # by hand for some values:
factor1 <- "TEs_rnaseq_class2"; factor2 <- "alltrue"; myacc <- 0.01
factor1 <- "srnas_unique_class2"; factor2 <- "alltrue"; myacc <- 0.01
factor1 <- "srnas_class2"; factor2 <- "alltrue"; myacc <- 0.01
factor1 <- "srnas_PenneC_class2"; factor2 <- "alltrue"; myacc <- 0.01


gglayers <- list(scale_fill_manual(values = TRIBE_colors),
                scale_shape_manual(values = F4_shapes),
                theme_bw(),
                theme(legend.position = "right",
                      axis.ticks = element_blank(),
                      axis.text.x = element_text(angle = 300, hjust = 0)))

myenrichment <- list()
myacc <- 0.1 # srnasrounding of y label
factor3     <-  "F4"
for(factor1 in factor1s) {
    #   factor1 <- "srnas_class2"
    #   factor1 <-  "meth_CHH_class_stringent"
  for (factor2 in factor2s) {
    message(factor2)
    #   factor2 <- factor2s[1]
    #   factor2 <- "alltrue"
    perc_df_all <- reduce(perc_df_all_lists[[factor1]][[factor2]], rbind)

    perc_df_filter <- perc_df_all %>%
      filter(!!sym(paste0(factor1))!= "inbetween") %>%
      filter(!!sym(paste0(factor1))!= "none") %>%
      filter(.[[2]] != "none")
    avg_de <- avg_de_filter[[factor1]]
    mymult <- n_distinct(perc_df_filter[,2])/2

    # calculate ratios in seperate df
    dat_text <- perc_df_filter %>%
      group_by(!!sym(paste0(factor2)), !!sym(paste0(factor1))) %>%
      summarise(
                #             n = n(),
                mean = mean(!!sym(paste0(factor2, "_freq")), na.rm = TRUE)) %>%
      filter(str_detect(!!sym(factor1), "^[+-].*")) %>%
      ungroup() %>%
      spread(!!sym(paste0(factor1)), mean)  %>%
      mutate(ratio = round(.[[3]] / .[[2]], 2))

     tmp <- list()
   # for (i in str_subset(levels(perc_df_filter[[1]]), "^[+-].*")) {
   for (i in str_subset((perc_df_filter[[1]]), "^[+-].*")) {
     # message(i)
     avg_med <- avg_de %>%
       filter(!!sym(paste0(factor1)) == i) %>%
       dplyr::select(median) %>%
       pull
       tmp[i] <- log2(dat_text[i]) - log2(avg_med)
     }
   tmp2 <- as.data.frame(sapply(tmp, function(x) x))
   tmp2$names <- as.character(dat_text[[1]])
   tmp3 <- pivot_longer(tmp2, cols = -names)

   gg4 <- ggplot(tmp3, aes(x = names, y = value, fill = name )) +
     geom_bar(stat="identity",position="dodge") +
     gglayers +
  labs(title = "",
  #    subtitle = "sub",
  #    caption="cap",
     x = "",
     y = "log(proportion / average proportion)")

   myenrichment[factor1] <- tmp3

    gg <- ggplot(perc_df_filter %>%
                 # filter our only DEG/DESL/DMR
                 filter(str_detect(!!sym(factor1), "^[+-].*")),
                 aes(x = !!sym(factor1), 
                     y = !!sym(paste0(factor2, "_freq")))) +
                facet_grid(vars(), vars(!!sym(factor2))) +
                scale_y_continuous(labels = scales::percent_format(acc = myacc, scale = 1),
                                   limits = c(0, NA)) +
                geom_hline(yintercept = avg_de$median[1], color = TRIBE_colors[5]) +
                geom_hline(yintercept = avg_de$median[2], color = TRIBE_colors[4]) +
                geom_boxplot(aes(fill = !!sym(factor1),alpha = 0.7)) +
                geom_point(aes(shape = !!sym(factor3))) +
                gglayers

    gg2 <- gg + geom_text(
                          data    = dat_text,
                          mapping = aes(x = -Inf, y = -Inf, label = ratio),
                          hjust   = "inward",
                          vjust   = -1
    ) 


    gg3 <- ggplot(perc_df_filter %>%
                 # filter our only DEG/DESL/DMR
                 filter(str_detect(!!sym(factor1), "^[+-].*")),
                 aes(x = !!sym(factor1), 
                     y = n)) +
                facet_grid(vars(), vars(!!sym(factor2))) +
                scale_y_continuous(
                                   limits = c(0, NA)) +
                    #   scale_y_continuous(labels = scales::unit_format(unit = "", scale = 1e-3, digits = 2),
                geom_boxplot(aes(fill = !!sym(factor1),alpha = 0.7)) +
                geom_point(aes(shape = !!sym(factor3))) +
                gglayers

                #                 ggsave(mytitle, gg, height = 6, width = mymult+2)
                mytitle <- paste0("stats_v2", "_allF4s_", factor1, "_", factor2)
                write_csv(perc_df_all, path = paste0(mytitle, ".csv"))
                ggall <- ggarrange(gg2, gg3)
                ggsave(paste0(mytitle, ".pdf"), ggall, height = 6, width = mymult*3+4, limitsize = F)
                ggsave(paste0(mytitle, "_enrichment.pdf"), gg4, height = 6, width = mymult+2, limitsize = F)


  }
}

# same for 2D
factor3     <-  "F4"
for(factor1 in factor1s) {
    #   factor1 <- "srnas_class2"
    #   factor1 <- "srnas_unique_class2"
  for (factor2a in factor2s) {
    for (factor2b in factor2s) {
      if (factor2b == factor2a) next
    #   factor2a <- factor2s[1]
    #   factor2b <- factor2s[3]
    message(factor2a,"_", factor2b)
    perc_df_all <- reduce(perc_df_all_lists_2d[[factor1]][[factor2a]][[factor2b]], rbind)

    perc_df_filter <- perc_df_all %>%
      filter(!!sym(paste0(factor1))!= "inbetween") %>%
      filter(!!sym(paste0(factor1))!= "none") %>%
      filter(.[[2]] != "none")
    avg_de <- avg_de_filter[[factor1]]
    mymult <- n_distinct(perc_df_filter[,2])
    mymult2 <- n_distinct(perc_df_filter[,3])

    gg <- ggplot(perc_df_filter %>%
                 # filter our only DEG/DESL/DMR
                 filter(str_detect(!!sym(factor1), "^[+-].*")),
                 aes(x = !!sym(factor1), 
                     y = !!sym(paste0(factor2a, "_freq")))) +
                                     facet_grid(cols = vars(!!sym(factor2a)), rows=vars(!!sym(factor2b))) +
                scale_y_continuous(labels = scales::percent_format(acc = 0.1, scale = 1)) +
                geom_hline(yintercept = avg_de$median[1], color = TRIBE_colors[5]) +
                geom_hline(yintercept = avg_de$median[2], color = TRIBE_colors[4]) +
                geom_boxplot(aes(fill = !!sym(factor1),alpha = 0.7)) +
                geom_point(aes(shape = !!sym(factor3))) +
                scale_shape_manual(values = F4_shapes) +
                scale_fill_manual(values = TRIBE_colors) +
                theme_bw() +
                theme(legend.position = "right",
                      axis.ticks = element_blank(),
                      axis.text.x = element_text(angle = 300, hjust = 0))

                mytitle <- paste0("stats_2d_", "_allF4s_", factor1, "_", factor2a, "_", factor2b, ".pdf")
                #                 ggsave(mytitle, gg, height = mymult2 + 2, width = mymult + 2, limitsize = F)
                #                 write_csv(perc_df_all, path = paste0(mytitle, ".csv"))

    gg3 <- ggplot(perc_df_filter %>%
                 # filter our only DEG/DESL/DMR
                 filter(str_detect(!!sym(factor1), "^[+-].*")),
                 aes(x = !!sym(factor1), 
                     y = n)) +
                    facet_grid(cols = vars(!!sym(factor2a)), rows=vars(!!sym(factor2b))) +
                    #                 scale_y_continuous(labels = scales::percent_format(acc = 0.1, scale = 1)) +
          scale_y_continuous(
                             breaks = scales::pretty_breaks(n = 8)) +
                geom_boxplot(aes(fill = !!sym(factor1),alpha = 0.7)) +
                geom_point(aes(shape = !!sym(factor3))) +
                scale_shape_manual(values = F4_shapes) +
                scale_fill_manual(values = TRIBE_colors) +
                theme_bw() +
                theme(legend.position = "right",
                      axis.ticks = element_blank(),
                      axis.text.x = element_text(angle = 300, hjust = 0))

                mytitle <- paste0("stats_2d_", "_allF4s_", factor1, "_", factor2a, "_", factor2b, ".pdf")

                ggall <- ggarrange(gg, gg3)
                ggsave(mytitle, ggall, height = 6, width = mymult*3+4, limitsize = F)

                #                 perc_df_filter2 <- perc_df_filter %>%
                #                 group_by(!!sym(factor2b), !!sym(factor1)) %>%
                #                   mutate(n_freq = n/sum(n)*300)
                # 
                #                 perc_df_filter2  %>%
                #                 group_by(!!sym(factor2a), !!sym(factor1)) %>%
                #                   summarize(nf= sum(n_freq), n = sum(n))
                # 
                #                 tmp <- perc_df_filter2  %>%
                #                   filter(TEs_Order_renamed == "TIR") %>%
                #                   filter(srnas_class2 == "-DESL")
                #             
                # 
                #                 group_by(!!sym(factor2a), !!sym(factor1)) %>%
                #                   summarize(nf= sum(n_freq), n = sum(n))

                #     gg4 <- ggplot(perc_df_filter2 %>%
                #                 filter our only DEG/DESL/DMR
                #                  filter(str_detect(!!sym(factor1), "^[+-].*")),
                #                  aes(x = !!sym(factor1), 
                #                      y = !!sym(paste0("n", "_freq")))) +
                #                                      facet_grid(cols = vars(!!sym(factor2a)), rows=vars(!!sym(factor2b))) +
                #                 scale_y_continuous(labels = scales::percent_format(acc = 0.1, scale = 1)) +
                #                 geom_hline(yintercept = avg_de$median[1], color = TRIBE_colors[5]) +
                #                 geom_hline(yintercept = avg_de$median[2], color = TRIBE_colors[4]) +
                #                 geom_boxplot(aes(fill = !!sym(factor1),alpha = 0.7)) +
                #                 geom_point(aes(shape = !!sym(factor3))) +
                #                 scale_shape_manual(values = F4_shapes) +
                #                 scale_fill_manual(values = TRIBE_colors) +
                #                 theme_bw() +
                #                 theme(legend.position = "right",
                #                       axis.ticks = element_blank(),
                #                       axis.text.x = element_text(angle = 300, hjust = 0))

  }
  }
}

# load(file="final_figures_seb_code2d.Rdata")
# load(file="final_figures_seb_code.Rdata")

#------------------------------------ {{{ targetting: extract sRNAs species from interessting bins:

  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter) %>%
    filter(srnas_unique_P4042_class2 == "+DESL")

  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter) %>%
    filter(srnas_P4042_class2 == "+DESL")

  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter) %>%
    filter(srnas_P4042_class2 == "-DESL")

  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter) %>%
    filter(srnas_P3611_class2 == "+DESL")


tmp=makeGRangesFromDataFrame(mydf_filtera, seqnames.field="#seqname")

sizefilter <- 22
p1 <- ScanBamParam( what=scanBamWhat())
path_map  <- file.path(path_srnas, "mapped_set3_4_SL30_Penn")
bam_F4s <- list.files(path_map, pattern = "*_KP....-1.*SolLyc.*.bam$")
bam <- "22_KP4042-1_SolLyc_zeroMM_bowtie.sorted.bam"
bam <- "19_KP3611-1_SolLyc_zeroMM_bowtie.sorted.bam"
bam <- "1_KM82C-1_SolLyc_zeroMM_bowtie.sorted.bam"

#      p2 <- ScanBamParam(which=mydf_filtera[1,], what=scanBamWhat())
readin <- readGAlignments(file.path(path_map, bam_F4s[1]), param=p1)
readin <- readin[width(readin) %in% sizefilter, ]
for(bam in bam_F4s[-1]){
  message(bam)
  tmp <- readGAlignments(file.path(path_map, bam), param=p1)
  readin <- c(readin, tmp[width(tmp) %in% sizefilter, ])
}

tmp2 <-      subsetByOverlaps(readin, tmp, ignore.strand=FALSE)
tmp3 <- ((DNAStringSet(tmp2@elementMetadata$seq)))
names(tmp3) <- tmp3

writeXStringSet(tmp3, "plusDESLP4042_multi.fa", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
writeXStringSet(tmp3, "minus_wt_DESLP4042_multi.fa", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
writeXStringSet(tmp3, "plusDESLP3611_multi.fa", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
# grep("Solyc11g0172",mydf_filtera$Genes_Slyc)

#}}}

#------------------------------------ importing sPARTA results {{{
srna_targets <- read_csv("/projects/TRIBE/srnas/sPARTA/P4042_top20sRNAs_plusDESL_predicted/All.targs.parsed.csv")
srna_targets <- read_csv("/projects/TRIBE/srnas/sPARTA/P3611_top20sRNAs_plusDESL_predicted/All.targs.parsed.csv")
srna_targets <- read_csv("/projects/TRIBE/srnas/sPARTA/P4042_top20sRNas_minusDESL_wt_predicted/All.targs.parsed.csv")
# how many target genes?
length(unique(srna_targets$Target))

myn <- 200

prefix <- "Genes_Slyc_rnaseq_"
load(file = file.path(path_comparison, "stats_general", paste0(prefix, "overlap_F4s.rdata"))) #list_deg, list_deg_overlap,list_deg_overlapNr,list_deg_overlapNrrandom, 
table(list_deg$P4042[,2])
table(list_deg$P4041[,2])

mydegs <- list_deg$P4042  %>%
  filter(rnaseq_P4042_class2 == "+DEG") %>%
  pull(1)

mydegs <- list_deg$P3611  %>%
  filter(rnaseq_P3611_class2 == "+DEG") %>%
  pull(1)

mydegs_minus <- list_deg$P4042  %>%
  filter(rnaseq_P4042_class2 == "-DEG") %>%
  pull(1)

mydegs_minus <- list_deg$P3611  %>%
  filter(rnaseq_P3611_class2 == "-DEG") %>%
  pull(1)

tmp3 <- replicate(myn,
mydegs_random <- list_deg$P4042  %>%
  filter(rnaseq_P4042_class2 == "non-DEG") %>%
  sample_n(length(mydegs_minus)) %>%
  pull(1)
  #   intersect(unique(tmp$Target), mydegs_random)
  #   intersect(unique(tmp$Target), as.vector(.$Genes_Slyc))
  #   print(.$Genes_Slyc)
  #   select(Genes_Slyc)
)


# intersect(substr(unique(srna_targets$Target),1,14), substr(mydegs, 1, 14))
# intersect(substr(unique(srna_targ-ts$Target),1,14), substr(mydegs_random, 1, 14))
sort(intersect(unique(srna_targets$Target), mydegs))
# [1] "Solyc07g039585.1.1" "Solyc07g041820.2.1" "Solyc07g038137.1.1" "Solyc12g062875.1.1" "Solyc08g048470.2.1" "Solyc08g080380.3.1" "Solyc01g104185.1.1"

sort(intersect(unique(srna_targets$Target), mydegs_minus))
# [1] "Solyc12g014500.2.1" "Solyc03g117600.3.1" "Solyc03g116630.3.1" "Solyc12g087840.1.1" "Solyc10g084890.2.1"
# for P3611
#  [1] "Solyc12g014500.2.1" "Solyc06g083720.2.1" "Solyc10g084670.2.1" "Solyc11g012860.2.1" "Solyc01g094690.3.1" "Solyc01g095360.2.1" "Solyc02g072240.3.1"
#  [8] "Solyc02g093980.3.1" "Solyc04g008820.3.1" "Solyc04g009050.3.1" "Solyc05g005180.3.1" "Solyc05g005760.3.1" "Solyc06g007440.3.1" "Solyc06g024380.1.1"
# [15] "Solyc06g070970.3.1"

myvec <- rep(NA,myn)
for(i in 1:myn) myvec[i] <- length((intersect(tmp3[,i], unique(srna_targets$Target))))
mean(myvec)
#---  }}}

#------------------- target analysis: working out profiles for slydf {{{
genes_sly_fasta <- "/data/public_data/tomato/ITAG3.2/ITAG3.2_CDS.fasta"
genes_fa <- readDNAStringSet(genes_sly_fasta, format="fasta")
# genes_fa2 <- genes_fa[lengths(genes_fa)>50 & alphabetFrequency(genes_fa)[,"N"]<1,]
genes_farev <- reverseComplement(genes_fa)

slydfmod <- slydf %>%
  #     filter(!!genotype == gtfilter) %>%
  #   select(ID, starts_with("srnas_P") & ends_with("class2"))
  #   dplyr::select(starts_with("srnas_P") & ends_with("class2")) %>%
  #   dplyr::filter(across((starts_with("srnas_P"))) == "DESL") 
  #     dplyr::filter(srnas_P1512_class == "DESL") 
  #     filter(srnas_P3611_class2 == "+DESL")
  mutate(across(matches("srnas_P.*class2"), 
                ~fct_recode(.x, "+" = "+DESL",
                            "-" = "-DESL",
                            "n" = "non-DESL",
                            "n" = "inbetween"))) %>%
  mutate(across(matches("rnaseq.*class2"), 
                ~fct_recode(.x, "+" = "+DEG",
                          "-" = "-DEG",
                          "n" = "non-DEG",
                          "n" = "inbetween")))

    myF4 ="P4042"
    infix = ""
    gtfilter = "00"
    myprefix <- "srnas_"
    suffix <- "_class2"
myparent <- "PenneC"
# myparent <- "M82C"

    # adding new columns:
    # n if hetero
    for (myF4 in F4s_srnas) {
    factor2  <- paste0(myprefix, myF4, infix, suffix)
    factor2b <- paste0(myprefix, myF4, infix, "_class3")
    factorg  <- paste0("genotype_", myF4)

      slydf <- slydf %>%
        # head(5000) %>%
        mutate(!!sym(factor2b) := case_when(
                  !!(sym(factorg)) != gtfilter ~ as.character(NA),
                  !!(sym(factor2)) == "-DESL" ~ "-",
                  !!(sym(factor2)) == "+DESL" ~ "+",
                  !!(sym(factor2)) == "non-DESL" ~ "n",
                  !!(sym(factor2)) == "inbetween" ~ as.character(NA) 
              ))
    }
    myprefix <- "rnaseq_"

    for (myF4 in F4s_rnaseq) {
    factor2  <- paste0(myprefix, myF4, infix, suffix)
    factor2b <- paste0(myprefix, myF4, infix, "_class3")
    factorg  <- paste0("genotype_", myF4)

      slydf <- slydf %>%
        # head(5000) %>%
        mutate(!!sym(factor2b) := case_when(
                  !!(sym(factorg)) != gtfilter ~ as.character(NA),
                  !!(sym(factor2)) == "-DEG" ~ "-",
                  !!(sym(factor2)) == "+DEG" ~ "+",
                  !!(sym(factor2)) == "non-DEG" ~ "n",
                  !!(sym(factor2)) == "inbetween" ~ as.character(NA) 
              ))

    }

        table2(slydf$srnas_P4042_class3)

        table2(slydf$TEs_rnaseq_P4042_class2,slydf$genotype_P4042)
        table2(slydf$TEs_rnaseq_P4042_class2,slydf$genotype_P4042)


slydfmod21 <- slydf %>%
  filter(srnas_2124ratio_class == "high_21bp")

cntsrna21 <- slydfmod21 %>%
  dplyr::count(across(matches("srnas_P.*class3")), sort = T)


#  head(cntsrna21, 20)
#   srnas_P1512_class2 srnas_P2561_class2 srnas_P2562_class2 srnas_P3611_class2 srnas_P3612_class2 srnas_P4041_class2 srnas_P4042_class2     n
# 1                  n                  n                  n                  n                  n                  n                  n 54681
# 2                  +                  +                  +                  +                  +                  +                  +  3252
# 3                  +                  +                  +                  n                  n                  n                  n  1445

  # rnaseq
df %>% mutate(across(c(x, starts_with("y")), mean, na.rm = TRUE))
nrow(distinct(tmp))
tmp2 <- head(tmp)
as.integer(tmp2[,2])

rowAny <- function(x) rowSums(x) > 1


cntrnaseq <- slydf %>% 
  distinct(Genes_Slyc, .keep_all= TRUE) %>%
  select(ID, Genes_Slyc, matches("^rnaseq_P.*class3$")) %>%
  filter(rowAny(across(matches("^rnaseq_P.*class3$"), ~ .x %in% c("+", "-")))) %>%
  dplyr::count(across(matches("^rnaseq_P.*class3$")), sort = T)
# Number 1-4, 7 and 20 look reasonable

write_csv(cntrnaseq, path = paste0("rnaseq_inheritance_profiles_v2.csv"))

  # rnaseq_P1512_class2 rnaseq_P2561_class2 rnaseq_P2562_class2 rnaseq_P3611_class2 rnaseq_P3612_class2 rnaseq_P4041_class2 rnaseq_P4042_class2     n
# 1                   n                   n                   n                   n                   n                   n                   n 25783
# 2                none                none                none                none                none                none                none  1711
# 3                   n                   n                   -                   n                   n                   n                   n   944
# 4                   n                   n                   +                   n                   n                   n                   n   799
# 5                   +                   n                   +                   n                   n                   n                   n   294


# P36* -DEG, +DESLs
tmp <- slydf %>%
  filter(across(matches("srnas_P36.*class2"), ~ (as.character(.x) == "-DESL" ))) %>%
  select((matches("srnas_P.*class2$")))

# srnas_prof_14 <- cntsrna21[14, ] %>% mutate(across(where(is.factor), as.character)) %>% unlist()
# rnaseq_prof_14 <- chartr("-+", "+-", srnas_prof_14)

# starting from rnaseq profile
## select a profile
rnaseq_prof_1 <- cntrnaseq[1, ] %>% mutate(across(where(is.factor), as.character)) %>% unlist()
rnaseq_prof_2 <- cntrnaseq[2, ] %>% mutate(across(where(is.factor), as.character)) %>% unlist()
rnaseq_prof_4 <- cntrnaseq[4, ] %>% mutate(across(where(is.factor), as.character)) %>% unlist()
# rnaseq_prof_7 <- cntrnaseq[7, ] %>% mutate(across(where(is.factor), as.character)) %>% unlist()
rnaseq_prof_19 <- cntrnaseq[19, ] %>% mutate(across(where(is.factor), as.character)) %>% unlist()
rnaseq_prof_20 <- cntrnaseq[20, ] %>% mutate(across(where(is.factor), as.character)) %>% unlist()
rnaseq_prof_23 <- cntrnaseq[23, ] %>% mutate(across(where(is.factor), as.character)) %>% unlist() # few if any
rnaseq_prof_26 <- cntrnaseq[26, ] %>% mutate(across(where(is.factor), as.character)) %>% unlist()# only 10 sRNA bins, no result
# starting from rnaseq profile
# sara: I like profiles 20 (4-), 24 (4+), and 27 (1not and 3 -).  Second choices: 8,9,12,15,17
head(cntrnaseq, 30)
i=19
rnaseq_prof <- cntrnaseq[i, ] %>% mutate(across(where(is.factor), as.character)) %>% unlist()
srnasprof <- chartr("-+", "+-", rnaseq_prof)


srnas_prof_fit <- slydf %>%
  filter(compareNA2(srnas_P1512_class3, srnasprof["rnaseq_P1512_class3"])) %>%
  filter(compareNA2(srnas_P2561_class3, srnasprof["rnaseq_P2561_class3"])) %>%
  filter(compareNA2(srnas_P2562_class3, srnasprof["rnaseq_P2562_class3"])) %>%
  filter(compareNA2(srnas_P3611_class3, srnasprof["rnaseq_P3611_class3"])) %>%
  filter(compareNA2(srnas_P3612_class3, srnasprof["rnaseq_P3612_class3"])) %>%
  filter(compareNA2(srnas_P4041_class3, srnasprof["rnaseq_P4041_class3"])) %>%
  filter(compareNA2(srnas_P4042_class3, srnasprof["rnaseq_P4042_class3"]))

srnas_prof_fit <- slydf %>%
  filter(compareNA(srnas_P1512_class3, srnasprof["rnaseq_P1512_class3"])) %>%
  filter(compareNA(srnas_P2561_class3, srnasprof["rnaseq_P2561_class3"])) %>%
  filter(compareNA(srnas_P2562_class3, srnasprof["rnaseq_P2562_class3"])) %>%
  filter(compareNA(srnas_P3611_class3, srnasprof["rnaseq_P3611_class3"])) %>%
  filter(compareNA(srnas_P3612_class3, srnasprof["rnaseq_P3612_class3"])) %>%
  filter(compareNA(srnas_P4041_class3, srnasprof["rnaseq_P4041_class3"])) %>%
  filter(compareNA(srnas_P4042_class3, srnasprof["rnaseq_P4042_class3"]))
# table(slydf$srnas_P1512_class3)

# srnas_prof_fit_nogenes <- srnas_prof_fit %>%
# rnaseq_prof <- rnaseq_prof_23 # few if any
#   filter(Genes_Slyc == "none")

# srnas_prof_fitting <- srnas_prof_fit_nogenes %>%
#   distinct(Genes_Slyc, .keep_all= TRUE) %>%
#   select(ID, Genes_Slyc, matches("^srnas_P.*class3$"))

table(srnas_prof_fit$srnas_2124ratio_class)
# high_24bp inbetween high_21bp
#       810      1567       809
table(srnas_prof_fit$dcl2_class2)
table(srnas_prof_fit$TEs_Order_renamed)
table(srnas_prof_fit$TEs_Superfamily)
table(srnas_prof_fit$annotation)

# rnaseq_prof_fit <- slydf %>%
#   filter(compareNA(rnaseq_P1512_class3, rnaseq_prof["rnaseq_P1512_class3"])) %>%
#   filter(compareNA(rnaseq_P2561_class3, rnaseq_prof["rnaseq_P2561_class3"])) %>%
#   filter(compareNA(rnaseq_P2562_class3, rnaseq_prof["rnaseq_P2562_class3"])) %>%
#   filter(compareNA(rnaseq_P3611_class3, rnaseq_prof["rnaseq_P3611_class3"])) %>%
#   filter(compareNA(rnaseq_P3612_class3, rnaseq_prof["rnaseq_P3612_class3"])) %>%
#   filter(compareNA(rnaseq_P4041_class3, rnaseq_prof["rnaseq_P4041_class3"])) %>%
#   filter(compareNA(rnaseq_P4042_class3, rnaseq_prof["rnaseq_P4042_class3"]))
rnaseq_prof_fit <- slydf %>%
  filter(compareNA2(rnaseq_P1512_class3, rnaseq_prof["rnaseq_P1512_class3"])) %>%
  filter(compareNA2(rnaseq_P2561_class3, rnaseq_prof["rnaseq_P2561_class3"])) %>%
  filter(compareNA2(rnaseq_P2562_class3, rnaseq_prof["rnaseq_P2562_class3"])) %>%
  filter(compareNA2(rnaseq_P3611_class3, rnaseq_prof["rnaseq_P3611_class3"])) %>%
  filter(compareNA2(rnaseq_P3612_class3, rnaseq_prof["rnaseq_P3612_class3"])) %>%
  filter(compareNA2(rnaseq_P4041_class3, rnaseq_prof["rnaseq_P4041_class3"])) %>%
  filter(compareNA2(rnaseq_P4042_class3, rnaseq_prof["rnaseq_P4042_class3"]))

# phasing bins:

srnas_prof_fit  <- slydf %>%
  filter(phasing_all22 != "none")

# extractin sRNAs expressed from fitting bins:
tmp5=makeGRangesFromDataFrame(srnas_prof_fit, seqnames.field="#seqname")
# now extracting sRNAs corresponding to those bins:
tmp3 <-      subsetByOverlaps(readin, tmp5, ignore.strand=FALSE) # 16k sRNAs
tmp4 <- ((DNAStringSet(tmp3@elementMetadata$seq)))
names(tmp4) <- tmp4
tmp6 <- sort(table(tmp4))
srnas_high_nr <- tmp6[as.integer(tmp6) > 100]
srnas_high_names <- names(srnas_high_nr)
length(srnas_high_nr) # 2k

# matches <- vmatchPDict(genes_farev,tmp4) # geht noch nicht (Biostrings 2.17.47)
# matches <- vcountPDict(genes_farev,tmp4)

matchesrev <- vcountPDict(PDict(DNAStringSet(srnas_high_names), tb.start = 9, tb.end = 12),
                       genes_farev, max.mismatch = 2)
matches <- vcountPDict(PDict(DNAStringSet(srnas_high_names), tb.start = 9, tb.end = 12),
                       genes_fa, max.mismatch = 2)

table(rowSums(matches))
table(colSums(matches))
table(rowSums(matchesrev))
table(colSums(matchesrev))

# inspect candidate genes/srnas pari
idxs <- which(colSums(matchesrev)>0)
genes_target <- names(genes_farev[idxs,])
genes_target_forward <- names(genes_farev[which(colSums(matches)>0),])


rnaseq_prof_fitting <- rnaseq_prof_fit %>%
  distinct(Genes_Slyc, .keep_all= TRUE) %>%
  select(ID, Genes_Slyc, matches("^rnaseq_P.*class3$"))

rnaseq_prof_target <- slydf %>%
  filter(substr(Genes_Slyc,1,14) %in% substr(genes_target,1,14)) %>%
  distinct(Genes_Slyc, .keep_all= TRUE) %>%
    select(ID, Genes_Slyc, matches("^rnaseq_P.*class3$"), rnaseq_dcl2_class2)
  # select(ID, Genes_Slyc, matches("^genotype.*"), matches("^rnaseq_P.*logFC$"))

rnaseq_prof_target_forward <- slydf %>%
  filter(substr(Genes_Slyc,1,14) %in% substr(genes_target_forward,1,14)) %>%
  distinct(Genes_Slyc, .keep_all= TRUE) %>%
    select(ID, Genes_Slyc, matches("^rnaseq_P.*class3$"), rnaseq_dcl2_class2)
  # select(ID, Genes_Slyc, matches("^genotype.*"), matches("^rnaseq_P.*logFC$"))

intersect(rnaseq_prof_target$Genes_Slyc, rnaseq_prof_fitting$Genes_Slyc)
print(i)
intersect(letters, "a")


gene="Solyc02g030137"
id=25 # print(genes_target) #and match
idxall <- idxs[id]
genes_farev[idxall,]
idxsrna <- which(matchesrev[,idxall] > 0)
mysrnas <- srnas_high_names[idxsrna]
srnas_high_nr[idxsrna]

matchesrev_sub <- vcountPDict(PDict(DNAStringSet(srnas_high_names[matchesrev[, idxall] > 0]), tb.start = 9, tb.end = 12),
                       genes_farev, max.mismatch = 3)
genes_target_sub <- names(genes_farev[which(colSums(matchesrev_sub)>0),])
readinsub <- readin[readin@elementMetadata$seq %in% mysrnas,]
tmp2 <- readinsub[!duplicated(start(readinsub)),]

tmp <- overlap_annot(tmp2, TEs[["Slyc"]], "TEs_Order", "Order", unique = FALSE)
tmp3 <-      subsetByOverlaps(sly_bins, tmp2, ignore.strand=FALSE)

as.data.frame(tmp3) %>%
  select(ID, matches("genotype"), TEs_Superfamily, matches("srnas_P.*class2$"))

# genes_target <- names(genes_fa[which(colSums(matches)>0),])
# genes_target <- names(genes_farev[which(colSums(matchesrev)>50),])
# DNAStringSet object of length 3:
#     width seq                                                                       names
# [1]  4629 TTAAGCTGAATTCCTCGATTTGCAAGAAAGAGTTG...GACCCCTCAATCTCAAGCTCGACAACGGTAGACAT Solyc01g008800.2....
# [2]   351 TTAGACATCCAACCCGTAATCGCAACGACCTAATT...ATGCACATCTCAAGGTCAGGTGCCGTTGAGCACAT Solyc11g021160.1....
# [3]   318 TCAAAGGGTAGAAGGGATATAGTGCATCAAGCTGT...TTCTTCGATTCTTTTTCCTGGCGCAGCTGAGCCAT Solyc11g044660.1....
# confirmed by psRNATarget that this targetting stratggy is a good approximation
srnas_high_names[which(rowSums(matchesrev)>0)]

# [38] "Solyc09g072690.1.1 Cation calcium exchanger (AHRD V3.3 *** A0A072U799_MEDTR)"
                   # ID         Genes_Slyc rnaseq_P1512_logFC rnaseq_P2561_logFC rnaseq_P2562_logFC rnaseq_P3611_logFC rnaseq_P3612_logFC rnaseq_P4041_logFC rnaseq_P4042_logFC
# 30 SL3.0ch09-65732201 Solyc09g072690.1.1              -5.37              -5.02              -7.15              -0.74              -3.66              -0.25              -0.58

[19] "Solyc04g026110.3.1 Disease resistance family protein (AHRD V3.3 *** B9GY55_POPTR)"
17 SL3.0ch04-18823801 Solyc04g026110.3.1                   n                   n                   -
17                   -                   -                   n                   n

vcountPDict(PDict(tmp6[1:3,], tb.start = 9, tb.end = 12), DNAStringSet(
            c("TAAGACTTAACTGGGTGACACA",
              "TAAGACTTAACTGGGTGACACAA",
              "CAAGACTTAACTGGGTGACACAA",
              "CAAGACTTAAATGGGTGACACAA",
              "TAAGACTTAACTGGGTGACACTA")),
            max.mismatch = 1)

matches <- vcountPDict(seqs_dict,submrnas)
matches[1185,4983]	# match
matchPattern(seqs1185,submrnas4983)

tmp %>% select(starts_with("srna")) %>% head

  select((matches("srnas_P.*class2$")))
#--- working out profiles for slydf }}}

#------------------- dcl2 DESL vs rnaseq {{{

factor1 <- "rnaseq_dcl2_class2_promoter"
factor1 <- "rnaseq_dcl2_class2"
factor2 <- "dcl2_class2"

mycounts <- slydf %>%
  dplyr::count(!!sym(factor1), !!sym(factor2)) %>%
  drop_na() %>%
  filter(!!sym(factor2) != "none") %>%
  group_by(!!sym(factor1)) %>%
  mutate(!!sym(paste0(factor1, "_freq")) := round(100 * n / sum(n),2)) %>%
  ungroup() %>%
  group_by(!!sym(factor2)) %>%
  mutate(!!sym(paste0(factor2,"_freq")) := round(100 * n / sum(n),2)) %>%
  ungroup()

    perc_df_filter <- mycounts %>%
      filter(!!sym(paste0(factor1))!= "inbetween") %>%
      filter(!!sym(paste0(factor1))!= "none") %>%
      filter(.[[2]] != "none")

    gg <- ggplot(mycounts %>%
                 # filter our only DEG/DESL/DMR
                 filter(str_detect(!!sym(factor1), "^[+-].*")),
                 aes(x = !!sym(factor1), 
                     y = !!sym(paste0(factor2, "_freq")))) +
                facet_grid(vars(), vars(!!sym(factor2))) +
                #                 scale_y_continuous(labels = scales::percent_format(acc = myacc, scale = 1),
                #                                    limits = c(0, NA)) +
                #                 geom_hline(yintercept = avg_de$median[1], color = TRIBE_colors[5]) +
                # geom_hline(yintercept = avg_de$median[2], color = TRIBE_colors[4]) +
                geom_boxplot(aes(fill = !!sym(factor1),alpha = 0.7)) +
                #                 geom_point(aes(shape = !!sym(factor3))) +
                #                 scale_shape_manual(values = F4_shapes) +
                #                 scale_fill_manual(values = TRIBE_colors) +
                theme_bw() +
                theme(legend.position = "right",
                      axis.ticks = element_blank(),
                      axis.text.x = element_text(angle = 300, hjust = 0))
                ggsave("dcl2.pdf", gg, height = 6, width = 9, limitsize = F)

    gg2 <- gg + geom_text(
                          data    = dat_text,
                          mapping = aes(x = -Inf, y = -Inf, label = ratio),
                          hjust   = "inward",
                          vjust   = -1
    ) 


    gg3 <- ggplot(perc_df_filter %>%
                 # filter our only DEG/DESL/DMR
                 filter(str_detect(!!sym(factor1), "^[+-].*")),
                 aes(x = !!sym(factor1), 
                     y = n)) +
                facet_grid(vars(), vars(!!sym(factor2))) +
                scale_y_continuous(
                                   limits = c(0, NA)) +
                    #   scale_y_continuous(labels = scales::unit_format(unit = "", scale = 1e-3, digits = 2),
                geom_boxplot(aes(fill = !!sym(factor1),alpha = 0.7)) +
                geom_point(aes(shape = !!sym(factor3))) +
                scale_shape_manual(values = F4_shapes) +
                scale_fill_manual(values = TRIBE_colors) +
                theme_bw() +
                theme(legend.position = "right",
                      axis.ticks = element_blank(),
                      axis.text.x = element_text(angle = 300, hjust = 0))

                #                 ggsave(mytitle, gg, height = 6, width = mymult+2)
                mytitle <- paste0("stats_v2", "_allF4s_", factor1, "_", factor2, ".pdf")
                write_csv(perc_df_all, path = paste0(mytitle, ".csv"))
                ggall <- ggarrange(gg2, gg3)
                ggsave(mytitle, ggall, height = 6, width = mymult*3+4, limitsize = F)

#--- dcl2 DESL vs rnaseq  }}}

#------------------- dcl2 targetting {{{

# which bins are DESL
dcldep <- sly_bins[slydf$dcl2_class2 == "-DESL", ]

dcldep$phasing_all22
tmp3 <-      subsetByOverlaps(readin, dcldep, ignore.strand=FALSE) # 16k sRNAs
tmp4 <- ((DNAStringSet(tmp3@elementMetadata$seq)))
names(tmp4) <- tmp4
tmp6 <- sort(table(tmp4))
# AAATACAATTTAGACGAAGCAA AAAATTTTATTCTAGATCCTAA
#                  24066                  37770
srnas_high_nr <- tmp6[as.integer(tmp6) > 1000]
srnas_high_names <- names(srnas_high_nr)
length(srnas_high_nr) # 140 sRNAs
 [1] "Solyc00g102400.3.1 NBS-LRR resistance protein-like protein"                                    "Solyc02g021480.3.1 Magnesium transporter CorA-like family protein (AHRD V3.3 *** AT1G29820.2)"
 [3] "Solyc05g043420.1.1 NBS-LRR resistance protein-like protein (AHRD V3.3 *** A1Y9R1_SOLLC)"       "Solyc05g054010.3.1 NBS-LRR resistance protein-like protein (AHRD V3.3 *** A1Y9R1_SOLLC)"
 [5] "Solyc06g008400.2.1 Mi1.5"                                                                      "Solyc06g008410.1.1 Disease resistance protein (AHRD V3.3 *-* Q9SBC3_SOLLC)"
 [7] "Solyc06g008450.3.1 Mi1.6"                                                                      "Solyc06g008480.2.1 Mi1.7"
 [9] "Solyc06g008770.2.1 CNL4"                                                                       "Solyc06g008800.2.1 CNL6"
[11] "Solyc06g048960.3.1 Dicer-like 2a"                                                              "Solyc06g065000.2.1 LOW QUALITY:Hero resistance protein (AHRD V3.3 *** Q8GSM1_SOLLC)"
[13] "Solyc11g008540.2.1 Dicer-like 2b"

Solyc00g102400.3....
seq <- "TATGCCACATCTAAAACACGTG"
                  7563

gene:Solyc02g021480.3 is up in DCL2 DEGS

!!!"Solyc06g048960.3.1 Dicer-like 2a"
1.5	4.58	422.361407007688	7.47242379151217E-94	0	inbetween	+	inbetween	191	222	207	692	672	679	gene:Solyc06g048960.3


readinsub <- readin[readin@elementMetadata$seq == seq,]
tmp2 <- readinsub[!duplicated(start(readinsub)),]

tmp <- overlap_annot(tmp2, TEs[["Slyc"]], "TEs_Order", "Order", unique = FALSE)
tmp3 <-      subsetByOverlaps(sly_bins, tmp2, ignore.strand=FALSE)
[1] Helitron_withPotentialHostGeneDomain Solyc06g008410.1.1        none                    Disease resistance protein (AHRD V3.3 *-* Q9SBC3_SOLLC) Solyc06g008410.1.1           none            00

# request
The first plot would be a genome-wide map of the high 22 EPRV DESLs.  
I would like to see whether they cluster at all (whether or not they do affects the design of the proposed experiments). The most useful format would be a set of plots showing where they are with all seven F4 high 22 EPRV DESLs on the same diagram (probably colour coded for each F4 and with separate + and – DESL plots).  Is there a chance you could get this one to me by 16thNovember or sooner. 
mydf_filtera %>%
  filter(srnas_P4042_class2 == "+DESL") %>%

slydf[,1] <- droplevels(slydf[,1])
    myprefix <- "srnas_"
for (myF4 in F4s_srnas) {

  #   myF4 <- F4s_rnaseq[1]
  message(myF4)
  genotype <- sym(paste0("genotype_", myF4))
  factor1  <- paste0(myprefix, myF4, infix, suffix)
  #   rnaseq_a  <- paste0(mydata, myF4, "_class")
  #save number of DEG/DET for each F4

  mydf_eprv <- slydf %>%
    filter(!!genotype == gtfilter) %>%
    filter(TEs_Superfamily == "EPRV") %>%
    distinct(TEs_name, .keep_all= TRUE)

  mydf_deslp <- slydf %>%
    filter(!!genotype == gtfilter) %>%
    filter(!!sym(factor1) == "+DESL") %>%
    filter(TEs_Superfamily == "EPRV") %>%
    distinct(TEs_name, .keep_all= TRUE)
    # filter(srnas_unique_P4042_class2 == "+DESL")
  mydf_deslm <- slydf %>%
    filter(!!genotype == gtfilter) %>%
    filter(!!sym(factor1) == "-DESL") %>%
    filter(TEs_Superfamily == "EPRV") %>%
    distinct(TEs_name, .keep_all= TRUE)
    # filter(srnas_unique_P4042_class2 == "+DESL")

  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter)

    # filter(srnas_2124ratio_class == "high_21bp")
  message(nrow(mydf_filtera))

  # Jan21: for degradome:D2G/+DEG
  # extract genes that are both up in F4 and dcl2 to see if they are  degradome targets

  # myF4 <- "P4041"
for (myF4 in F4s_rnaseq) {
  genotype <- sym(paste0("genotype_", myF4))
  # genotype_P2562
  tmp <- slydf %>%
    filter(!!genotype == gtfilter) %>%
    filter(!!sym(paste0("rnaseq_",myF4,"_class2")) == "+DEG") %>%
    filter(rnaseq_dcl2_class2 == "+D2G") %>%
    # distinct(Genes_Slyc, .keep_all= TRUE)
    distinct(Genes_Slyc) %>%
    mutate(Genes_Slyc = str_replace(Genes_Slyc, "\\..*",""))
  write_csv(tmp, path = paste0("UPdeg",myF4,"upD2G_seb", ".csv"))
}
    
#--- dcl2 targetting  }}}

#------------------- attempt of genomewide coverage plots {{{
ggplot(mydf_deslp, aes(aesthetics))

grgt=makeGRangesFromDataFrame(mydf_filtera, seqnames.field="#seqname")
seqlengths(grgt) <- width(grtomato)
grgt <- IRanges::reduce(grgt)

grdeslp=makeGRangesFromDataFrame(mydf_deslp, seqnames.field="#seqname")
seqlengths(grdeslp) <- width(grtomato)

grdeslm=makeGRangesFromDataFrame(mydf_deslm, seqnames.field="#seqname")
seqlengths(grdeslm) <- width(grtomato)

greprv=makeGRangesFromDataFrame(mydf_eprv, seqnames.field="#seqname")
seqlengths(greprv) <- width(grtomato)

seqlengths(tmp2) <- seqlengths(tmp2)[1:13]
seqlengths(tmp2) <- seqlengths(tmp2)[1:13]
tmp2@seqnames <- droplevels(tmp2@seqnames)
tmp2@seqnames <- grtomato@seqnames
table(mydf_filtera[,1])

      #seqname    start      end width strand                 ID       genome type chromatin_state
3634 SL3.0ch12 60413601 60413800   200      * SL3.0ch12-60413601 SL3.0-Genome  200            none

  library(ggbio)

tomato <- Rle(levels(droplevels(seqnames(sly_bins))), (rep(1,13)))

tomato_fa <- readDNAStringSet(genome_sly_path, format="fasta")

grtomato <- GRanges(seqnames=tomato,
IRanges(1, width=(seqlengths(sly_bins)[1:13])),
stand="*")
seqlengths(grtomato) <- width(grtomato)
grtomato@seqnames
tmp2@seqnames

p <- ggbio() +
circle(grgt, geom = "rect", color = "steelblue") +
circle(grdeslp, geom = "heatmap", color = "green", size = 0.5, scale.n = 5000) +
# circle(grdeslp, geom = "rect",  trackWidth = 2, size = 0.5, scale.n = 5000) +
circle(grtomato, geom = "ideo", fill = "gray70") +
circle(grtomato, geom = "scale", size = 2) +
circle(grtomato, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
# circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
pdf("tmp.pdf")
print(p)
# autoplot(grdeslp, layout = "karyogram")
dev.off()


## average over gr object
sly_bins50k <- subset(bins[["50k"]]$gr, genome == "SL3.0-Genome")
sly_bins50k$deslp <- countOverlaps(sly_bins50k, grdeslp)
sly_bins50k$deslm <- -countOverlaps(sly_bins50k, grdeslm)
sly_bins50k$eprv <- countOverlaps(sly_bins50k, greprv)
tmp <- as.data.frame(sly_bins50k)

gg <- ggplot(tmp, aes(x = start, y = deslp, color = "blue")) +
  geom_line() +
  geom_line(aes(y = eprv) ,color = "grey") +
  geom_line(aes(y = deslm) ,color = "steelblue") +
  facet_grid(seqnames~.)
ggsave("genomewide.pdf", gg, height = 8, width = 9, limitsize = F)
#--- attempt of genomewide coverage plots }}}


#------------------- comparing DET and DESL for good candiates {{{

myF4 ="P4042"
infix = ""
suffix <- "_class2"
# gtfilter = "11" # for penn homozygous
# gtfilter = "00" # for lyc
# myprefix <- "srnas_unique_"
myprefix <- "srnas_"
factor1  <- paste0(myprefix, myF4, infix, suffix)
myprefix <- "TEs_rnaseq_"
factor2  <- paste0(myprefix, myF4, infix, suffix)
# [1] "TEs_rnaseq_P4042_class2"
factor1

slydf %>%
  count(!!sym(factor2), TEs_Order_renamed, !!sym(factor1)) %>%
  pivot_wider(names_from = !!sym(factor2), values_from = n) %>%
  arrange(desc(!!sym(factor1))) %>%
  as.data.frame()


tmp <- with(slydf, table(TEs_rnaseq_P4042_class2, TEs_Order_renamed))
round(proportions(tmp,2)*100,1)
#                        TEs_Order_renamed
# TEs_rnaseq_P4042_class2   TIR   LTR  none  LINE Pararetrovirus Helitron Other_TEs
#               -DET        0.1   0.1   0.0   0.3            0.1      0.5       0.1
#               +DET        0.2   0.2   0.0   0.2            0.5      0.7       0.1
#               inbetween   1.0   1.0   0.0   1.0            2.4      1.5       0.5
#               non-DET    18.3  19.9   0.0  18.7           16.0     23.7      15.9
#               none       80.5  78.9 100.0  79.7           80.9     73.6      83.5
tmp <- with(slydf, table(srnas_P4042_class2, TEs_Order_renamed))
round(proportions(tmp,2)*100,1)
#                   TEs_Order_renamed
# srnas_P4042_class2  TIR  LTR none LINE Pararetrovirus Helitron Other_TEs
#          -DESL      3.9  0.7  1.8  1.4            3.1      1.1       2.7
#          +DESL      1.3  0.9  0.6  0.5           13.6      2.4       0.7
#          inbetween 37.8 28.4 25.1 30.7           40.5     32.5      32.0
#          non-DESL  57.0 70.1 72.6 67.3           42.8     64.0      64.6
tmp <- with(slydf, table(srnas_P4042_class2, TEs_Order_renamed))
round(proportions(tmp,2)*100,1)


tmp <- with(slydf, table(srnas_P4042_class2, TEs_Order_renamed))

#                        TEs_Order_renamed
# TEs_rnaseq_P4042_class2     TIR     LTR    none    LINE Pararetrovirus Helitron Other_TEs
#               -DET          112    1483       0     153             41       41       519
#               +DET          374    2677       0     134            141       66       511
#               inbetween    1757   14920       0     568            695      134      2975
#               non-DET     32355  308570       0   10285           4551     2121     97407
#               none       142688 1222997 1348882   43816          23025     6576    512600
# within LTRs mostly Copia
with(slydf, table(srnas_P4042_class2, TEs_rnaseq_P4042_class2, TEs_Order_renamed))
# , , TEs_Order_renamed = Pararetrovirus
                  # TEs_rnaseq_P4042_class2
# srnas_P4042_class2   -DET   +DET inbetween non-DET   none
         # -DESL          0      0        14      86    770
         # +DESL         11     49       137     504   3178
         # inbetween     21     50       276    1754   9434
         # non-DESL       9     42       268    2207   9643
# , , TEs_Order_renamed = LTR
                  # TEs_rnaseq_P4042_class2
# srnas_P4042_class2   -DET   +DET inbetween non-DET   none
         # -DESL         41     12        75    1052   9350
         # +DESL         25    126       662    2466  10221
         # inbetween    571   1026      5243   84132 349110
         # non-DESL     846   1513      8940  220920 854316
with(slydf, table(srnas_unique_P3612_class2, TEs_rnaseq_P3612_class2, TEs_Order_renamed))
with(slydf, table(srnas_P4042_class2, TEs_Order_renamed))
with(slydf, table(factor1, factor2, TEs_Order_renamed))
#--- comparing DET and DESL for good candiates }}}

#------------------- heatmaps to track sRNAs over generations {{{
library("pheatmap")

load(file=file.path(path_srnas, "counts_cpm_alllineages_srnas_tribe.rdata")) # cnt_filt

lineages <- unique(str_replace(colnames(cnt_filt), "-.*", ""))
colnames(cnt_filt) <- str_replace(colnames(cnt_filt), "-", "x")

mylin <- lineages[1]

for (mylin in lineages) {
  message(mylin)
  try(
  cnt_filt <- cnt_filt %>%
    # mutate(test = (srnas_P4042x1_cpm))
    mutate(!!sym(paste0(mylin, "_avg_cpm")) := (!!sym(paste0(mylin, "x1_cpm"))+!!sym(paste0(mylin, "x2_cpm")))/2)
  )
}

tmp <- cnt_filt %>%
  dplyr::select(!contains("x"))

slydf <- cbind(slydf, tmp)
# save slydf v38

filt <- slydf %>%
  slice_max(srnas_P4042_22bp_count, n = 40) %>%
  dplyr::select(matches("srnas_P4.*_cpm") | matches("srnas_M.*_cpm") ) %>%
  
myF4="P4042"
genotype <- sym(paste0("genotype_", myF4))
infix=""
factor1  <- paste0(myprefix, myF4, infix, suffix)

myF4="P3611"
pat="srnas_P3.*_cpm"
genotype <- sym(paste0("genotype_", myF4))

myF4="P2561"
pat="srnas_P2.*_cpm"
genotype <- sym(paste0("genotype_", myF4))

myF4="P1512"
pat="srnas_P1.*_cpm"
genotype <- sym(paste0("genotype_", myF4))

myprefix <- "srnas_unique_"
myprefix <- "srnas_"
infix=""
suffix <- "_class2"
factor1  <- paste0(myprefix, myF4, infix, suffix)
suffix <- "_FDR"
factor_fdr  <- paste0(myprefix, myF4, infix, suffix)
suffix <- "_logFC"
factor_lfc  <- paste0(myprefix, myF4, infix, suffix)

sign="+DESL"
sign="-DESL"
  
# filt <- slydf2 %>%
#     filter(!!genotype == gtfilter) %>%
#     filter(!!sym(factor1) == "-DESL") %>%
#     filter(TEs_Superfamily == "EPRV")
#     # slice_max(srnas_P4042_22bp_count, n = 90) %>%

filt <- slydf %>%
    filter(!!genotype == gtfilter) %>%
    # filter(!!sym(factor1) == sign) %>%
    filter(!!sym(factor_fdr) < 1e-10) %>%
    filter(abs(!!sym(factor_lfc)) > 9) %>%
    # filter(TEs_Superfamily == "EPRV") %>%
    dplyr::select(matches("srnas_M.*_cpm") | matches(pat) | TEs_Order_renamed | ID | srnas_dcl2_class2 | annotation, !!sym(factor1)) %>%
    arrange(desc(TEs_Order_renamed)) %>%
    arrange(desc(!!sym(factor1)))


# seperate between data and annotation:
filtnum <- filt %>%
  dplyr::select(where(is.numeric))
filtannot <- filt %>%
  dplyr::select(-c(where(is.numeric) | ID))
rownames(filtannot) <- rownames(filtnum) <- filt$ID

# mat <- filt %>%
#   dplyr::select(matches("srnas_P4.*_cpm") | matches("srnas_M.*_cpm") )

# df <- filt %>%
#   dplyr::select(matches("srnas_P4042_class2") | matches("annotation") )
pdf(file=paste0("heatmap_lineages", myF4, ".pdf"),onefile=F, width = 11, height = 29) # or other device
pheatmap(filtnum, cluster_rows = F, cluster_cols = FALSE, scale="row", annotation_row = filtannot, cutree_rows = 2)
dev.off()

# for +DESL
if (sign == "+DESL") {
filtn <- filtnum %>%
    mutate(across(where(is.numeric), ~ .x / srnas_P4042_avg_cpm))
} else {
filtn <- filtnum %>%
    mutate(across(where(is.numeric), ~ .x / srnas_M82C_avg_cpm))
}


pdf(file=paste0("heatmap_lineages_", sign , myF4, ".pdf"),onefile=F, width = 11, height = 29) # or other device
pheatmap(filtn, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = filtannot)
dev.off()


gg <- pheatmap(filt, annotation_col=df)

pdf(file="heatmap_lineages.pdf",onefile=F, width = 16, height = 8) # or other device
pheatmap(filtn, cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(filtn, cluster_rows = F, cluster_cols = FALSE, annotation_row = filtannot, cutree_rows = 2)
dev.off()
# gg <- pheatmap(mat, annotation_col=df)

# pdf(file="heatmap_lineages.pdf",onefile=F, width = 16, height = 8) # or other device

# gg <- pheatmap(mat,
# scale="row",
# cluster_rows = FALSE, 
# cluster_cols = FALSE)
# dev.off()

#--- heatmaps to track sRNAs over generations }}}

#------------------- DMR vs EPRV {{{
setwd(file.path(path_srnas, "srnas_vs_meth_eprv"))
# May21: revisiting DMRs to ask the question about if
# EPRVs are activated by methylation changes
# approach: check if there are DMRs overlapping with EPRVs worth checking out
# reload v34 (which includes meth data)
load(file = file.path(path_comparison, paste0(mygenome, "_slydf_plus_statsv34.rdata")))

    myprefix <- "meth_"
    suffix <- "_class_stringent"
    infix <- "CHH"
    factor1  <- paste0(myprefix, myF4, "_", infix, suffix)

    myprefix <- "srnas_"
    suffix <- "_class2"
    infix = ""
    factor2  <- paste0(myprefix, myF4, infix, suffix)

factor1s <- c("srnas_class2", "srnas_unique_class2", "TEs_rnaseq_class2", "rnaseq_class2", "meth_CHH_class_stringent", "meth_CpG_class_stringent", "meth_CHG_class_stringent") 
    #   factor1 <-  "meth_CHH_class_stringent"

  for (infix in contexts) {
    # infix = "CpG"
    # factor1  <- paste0(myprefix, myF4, "_", infix, suffix)
  }

for (myF4 in F4s_srnas) {
  #   myF4 <- F4s_srnas[7]
  message(myF4)
  genotype <- sym(paste0("genotype_", myF4))
    myprefix <- "srnas_"
    suffix <- "_class2"
    infix = ""
  factor1  <- paste0(myprefix, myF4, infix, suffix)
  #   rnaseq_a  <- paste0(mydata, myF4, "_class")
  #save number of DEG/DET for each F4
  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter)
  message(nrow(mydf_filtera))

  try(
  for (infix in contexts) {
    myprefix <- "meth_"
    suffix <- "_class_stringent"
    factor1  <- paste0(myprefix, myF4, "_", infix, suffix)

    tmp <- mydf_filtera %>%
      # filter(srnas_P4042_class2 == "+DESL") %>%
      # dplyr::count(!!sym(factor1), !!sym(factor2)) %>%
      # drop_na() %>%
      filter(!!sym(factor2) %in% c("-DESL", "+DESL")) %>%
      filter(!!sym(factor1) %in% c("-DMR", "+DMR")) %>%
      filter(TEs_Order_renamed == "Pararetrovirus") %>%
      # group_by(!!sym(factor1)) %>%
      # mutate(!!sym(paste0(factor1, "_freq")) := round(100 * n / sum(n),2)) %>%
      dplyr::rename( "seqnames" = "#seqname")
    write_tsv(tmp, file = paste0("EPRV_DESL_DMR_",infix, "_",myF4,"_list_bins", ".bed"))

    tmp2 <- mydf_filtera  %>%
      # filter(!!sym(factor2) %in% c("-DESL", "+DESL")) %>%
      filter(!!sym(factor1) %in% c("-DMR", "+DMR")) %>%
      # dplyr::count(!!sym(factor1), !!sym(factor2),TEs_Order_renamed)
      dplyr::count(!!sym(factor1), TEs_Order_renamed)
    write_tsv(tmp2, file = paste0("DESL_DMR_",infix,"_",myF4,"_stats_bins", ".bed"))
  }
  )
}
# not many DMR/DESL/EPRV bins were found above (only 1-2 per F4)
# new strategy:
# - ranking EPRVs that are "unique", i.e. having a high unique cpm

srnas_unique_P4042_logCPM

# setwd(file.path(path_srnas, "srnas_vs_meth_eprv")) 
for (myF4 in F4s_srnas) {
  genotype <- sym(paste0("genotype_", myF4))
  myprefix <- "srnas_unique_"
  suffix <- "_logCPM"
  infix = ""
  factor1  <- paste0(myprefix, myF4, infix, suffix)

  myprefix <- "srnas_"
  suffix <- "_class"
  suffix <- "_class2"
  factor2  <- paste0(myprefix, myF4, infix, suffix)
  #   rnaseq_a  <- paste0(mydata, myF4, "_class")
  #save number of DEG/DET for each F4
  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter)
  message(nrow(mydf_filtera))

  print(factor1)
  print(factor2)
  # [1] "srnas_unique_P4042_logCPM"
  rank_list <- mydf_filtera  %>%
    filter(TEs_Order_renamed == "Pararetrovirus") %>%
    filter(!!sym(factor2) %in% c("-DESL", "+DESL")) %>%
    filter(!!sym(factor1) > 1) %>%
    add_count(TEs_name) %>%
    dplyr::rename( "number_bins_per_TE" = "n") %>%
    distinct(TEs_name, .keep_all = T) %>%
    dplyr::rename( "seqnames" = "#seqname") %>%
    arrange(!!sym(factor1))

  write_tsv(rank_list, file = paste0("EPRV_DESL_",myF4,"_ranking", ".bed"))
}

    # arrange(srnas_unique_P4042_logCPM)
    # pull(srnas_P4042_logCPM)

# adding TEs promotors
load(file = file.path(path_comparison, paste0("R-objects-bins-SL30-Penn_v7_annot.rdata")))
Promoters_TEs <- promoters(TEs[[ "Slyc" ]], upstream=500, downstream=50, use.names=TRUE)
TEs[["Slyc"]]$length <- width(TEs[["Slyc"]])
export.gff3(Promoters_TEs, con = file.path(paste0("Promoters_TEs_Slyc300", ".gff")))
sly_bins <- subset(bins[["200"]]$gr, genome == "SL3.0-Genome")
# tmp <- overlap_annot(sly_bins, TEs[["Slyc"]], "Promoters_TEs_Order", "Order", unique = FALSE)
sly_bins <- overlap_annot(sly_bins, TEs[["Slyc"]], "TEs_length", "length", unique = FALSE)
sly_bins <- overlap_annot(sly_bins, Promoters_TEs, "Promoters_TEs_Order", "Order", unique = FALSE)
sly_bins <- overlap_annot(sly_bins, Promoters_TEs, "Promoters_TEs_Name", "Name", unique = FALSE)
sly_bins <- overlap_annot(sly_bins, Promoters_TEs, "Promoters_TEs_Class", "Class", unique = FALSE)
sly_bins <- overlap_annot(sly_bins, Promoters_TEs, "Promoters_TEs_Superfamily", "Superfamily", unique = FALSE)
# add to big table
slydf$TEs_length <- sly_bins$TEs_length
slydf$Promoters_TEs_Order <- sly_bins$Promoters_TEs_Order
slydf$Promoters_TEs_Name <- sly_bins$Promoters_TEs_Name
slydf$Promoters_TEs_Class <- sly_bins$Promoters_TEs_Class
slydf$Promoters_TEs_Superfamily <- sly_bins$Promoters_TEs_Superfamily
save(slydf, file = file.path(path_comparison, paste0(mygenome, "_slydf_plus_statsv39.rdata")))


for (myF4 in F4s_srnas) {
    myprefix <- "meth_"
    suffix <- "_class_stringent"
    infix <- "CpG"
    factor1  <- paste0(myprefix, myF4, "_", infix, suffix)
# do EPRV DESLs have more DMRs?
  try(
    tmp <- slydf  %>%
      filter(TEs_Order_renamed == "Pararetrovirus") %>%
      # filter(!!sym(factor1) %in% c("-DMR", "+DMR")) %>%
      count(!!sym(factor1))
  )
  print(tmp)
}
# table is in wiki

  for (infix in contexts) {
    print(infix)
for (myF4 in F4s_meth2) {
    myprefix <- "meth_"
    suffix <- "_class_stringent"
    # infix <- "CpG"
    factor1  <- paste0(myprefix, myF4, "_", infix, suffix)
  myprefix <- "srnas_"
  suffix <- "_class2"
  factor2  <- paste0(myprefix, myF4, "", suffix)
    tmp <- slydf  %>%
      filter(TEs_Order_renamed == "Pararetrovirus") %>%
      filter(!!sym(factor2) %in% c("-DESL", "+DESL")) %>%
      filter(!!sym(factor1) %in% c("-DMR", "+DMR")) %>%
      # filter(!!sym(factor1) != "inbetween") %>%
      # filter(!!sym(factor2) != "inbetween") %>%
      count(!!sym(factor1), !!sym(factor2))
    print(tmp)
}
  }

  # most meth_P1512_CHG_class_stringent
factor1="meth_P2561_CHG_class_stringent"
factor2="srnas_P2561_class2"
tmp <- slydf  %>%
  filter(TEs_Order_renamed == "Pararetrovirus") %>%
  filter(!!sym(factor2) %in% c("-DESL", "+DESL")) %>%
  filter(!!sym(factor1) %in% c("-DMR", "+DMR")) %>%
  # filter(!!sym(factor1) != "inbetween") %>%
  # filter(!!sym(factor2) != "inbetween") %>%
  select(TEs_name, matches("^meth_P.*stringent$"))

write_csv(tmp, path = paste0(factor1, factor2, "_table", "_", ".csv"))

#--- DMR vs EPRV }}}
tt=sly_bins[1,]$TEs_name
idx=which(TEs[["Slyc"]]$Name==tt)
TEs[["Slyc"]][idx,]
