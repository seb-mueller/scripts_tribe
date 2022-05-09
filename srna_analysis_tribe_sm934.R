#------------------------------------ init {{{
library(readr)
library(edgeR)
library(dplyr)
path_base <-  "/projects/TRIBE"
path_base <-  "~/workspace/tribe/"
source(file.path(path_base, "scripts/project_settings.R"))
load(file = file.path(path_comparison, paste0(mygenome, "_slydf_plus_statsv39.rdata")))

colnames(slydf)

# slydf <- slydf  %>%
  # dplyr::select(!contains("meth")) 
# load(file = file.path(path_comparison, paste0("R-objects-bins-SL30-Penn_v6_annot.rdata")))

path_comparison <- file.path(path_base, "comparisons_using_big_table/")
mygenome  <- "SolLyc"
# load(file = file.path(path_comparison, paste0(mygenome, "_bins_lyc.rdata")))

prefix <- "sRNA_bins200_"
path_srnas  <- file.path(path_base, "srnas")
path_map  <- file.path(path_srnas, "mapped_set3_4_SL30_Penn")
# path_map <- "/projects/TRIBE/srnas/mapped_set3_4_SL30_Penn"
path_map_unique <- file.path(path_map, "mapped_unique_merged")
path_map_dcl2 <- file.path(path_srnas, "dcl2_zhengming")
# [1] "/projects/TRIBE/srnas/mapped_set3_4_SL30_Penn/mapped_set3_4_SL30_Penn/mapped_unique_merged"


meta <- read_csv(file.path(path_srnas, "sRNA_metadata_sara.csv"))
meta_dcl2 <- read_csv(file.path(path_srnas, "sRNA_metadata_dcl2.csv"))
meta_unique <- read_csv(file.path(path_srnas, "sRNA_metadata_sara_unique.csv"))
# mybins
setwd(file.path(path_srnas, "analysis_bins/"))
setwd(file.path(path_srnas))
# setwd(path_map)


meta_sub <- meta %>%
  filter(include == "yes") %>%
  # filter(plant == "P3611") %>%
  filter(mapped == mygenome)

# all are mapped == "merged"
meta_unique_sub <- meta_unique %>%
  filter(include == "yes")

# fls <- list.files(path_map, pattern=".bam$", full.names =F)
# fls_full <- list.files(path_map, pattern=".bam$", full.names =T)
# meta2 <- meta %>%
#   mutate(id = str_split(File,"_") %>% map_chr(.,1)) %>%
#   mutate(lib = str_split(File,"_") %>% map_chr(.,2)) %>%
#   mutate(plant = str_sub(str_split(lib,"-") %>% map_chr(.,1), 2)) %>%
#   mutate(rep = str_split(lib,"-") %>% map_chr(.,2)) %>%
#   mutate(mapped = str_split(File,"_") %>% map_chr(.,3))
# write_csv(meta2, "sRNA_metadata.csv")
fls <- meta_sub$File
fls_unique <- meta_unique_sub$File

# }}}

#------------------------------------ read in sRNA counts per size (cnts, adds to slydf){{{
#reading in sRNAs raw data and overlap with bins
  # counts <- gwcov2(mybins, meta[1:2,], path = path)
# mybins_sub is most likely sly_bins!
mycounts <- gwcov2(sly_bins, as.data.frame(meta_sub), path = path_map)

mycounts_dcl2 <- gwcov2(sly_bins, as.data.frame(meta_dcl2), path = path_map_dcl2)$wide
mycounts20dcl2 <- gwcov2(sly_bins, as.data.frame(meta_dcl2), path = path_map_dcl2, sizefilter = c(20))$wide
mycounts21dcl2 <- gwcov2(sly_bins, as.data.frame(meta_dcl2), path = path_map_dcl2, sizefilter = c(21))$wide
mycounts22dcl2 <- gwcov2(sly_bins, as.data.frame(meta_dcl2), path = path_map_dcl2, sizefilter = c(22))$wide
mycounts23dcl2 <- gwcov2(sly_bins, as.data.frame(meta_dcl2), path = path_map_dcl2, sizefilter = c(23))$wide
mycounts24dcl2 <- gwcov2(sly_bins, as.data.frame(meta_dcl2), path = path_map_dcl2, sizefilter = c(24))$wide
mycounts25dcl2 <- gwcov2(sly_bins, as.data.frame(meta_dcl2), path = path_map_dcl2, sizefilter = c(25))$wide
save(mycounts20dcl2,mycounts23dcl2,mycounts25dcl2,mycounts21dcl2,mycounts22dcl2,mycounts24dcl2,mycounts_dcl2,file = paste0(prefix, mygenome, "_counts_dcl2.rdata"))
# load(file = paste0(prefix, mygenome, "_counts_dcl2.rdata"))

mycounts_unique <- gwcov2(sly_bins, as.data.frame(meta_unique), path = path_map_unique)
# mycounts21 <- gwcov2(sly_bins, as.data.frame(meta_sub), path = path_map, sizefilter = c(21,22))
mycounts20only <- gwcov2(sly_bins, as.data.frame(meta_sub), path = path_map, sizefilter = 20)$wide
mycounts21only <- gwcov2(sly_bins, as.data.frame(meta_sub), path = path_map, sizefilter = c(21))$wide
mycounts22only <- gwcov2(sly_bins, as.data.frame(meta_sub), path = path_map, sizefilter = c(22))$wide
mycounts23only <- gwcov2(sly_bins, as.data.frame(meta_sub), path = path_map, sizefilter = c(23))$wide
mycounts24only <- gwcov2(sly_bins, as.data.frame(meta_sub), path = path_map, sizefilter = 24)$wide
mycounts25only <- gwcov2(sly_bins, as.data.frame(meta_sub), path = path_map, sizefilter = 25)$wide
save(mycounts20only,mycounts25only,  mycounts21only, mycounts22only, mycounts23only, mycounts24only, file = paste0(prefix, mygenome, "_counts_misc_filter.rdata"))
# rm(mycounts20only,mycounts25only,  mycounts21only, mycounts22only, mycounts23only, mycounts24only)
load(file = paste0(prefix, mygenome, "_counts_misc_filter.rdata"))
# save(mycounts, mycounts21, mycounts24only, file = paste0(prefix, mygenome, "_counts_misc.rdata"))
load(file = file.path(path_srnas, paste0(prefix, mygenome, "_counts.rdata")))
cnts <- mycounts$wide
# cnts_dcl2 <- mycounts_dcl2$wide
cnts_unique <- mycounts_unique$wide
cnts21 <- rowSums(mycounts21$wide[, meta_sub$generation == "F4"])
cnts24 <- rowSums(mycounts24only$wide[, meta_sub$generation == "F4"])
# cnts21only <- rowSums(mycounts21only$wide[, meta_sub$generation == "F4"])
# cnts22only <- rowSums(mycounts22only$wide[, meta_sub$generation == "F4"])
# cnts23only <- rowSums(mycounts23only$wide[, meta_sub$generation == "F4"])
# cnts24only <- rowSums(mycounts24only$wide[, meta_sub$generation == "F4"])
# cnts21only <- rowSums(mycounts21only$wide[, which(meta_sub$plant == "M82C")])
# cnts22only <- rowSums(mycounts22only$wide[, which(meta_sub$plant == "M82C")])
# cnts23only <- rowSums(mycounts23only$wide[, which(meta_sub$plant == "M82C")])
# cnts24only <- rowSums(mycounts24only$wide[, which(meta_sub$plant == "M82C")])
# cnts21only <- rowSums(mycounts21only$wide[, 1:2])
# cnts22only <- rowSums(mycounts22only$wide[, 1:2])
# cnts23only <- rowSums(mycounts23only$wide[, 1:2])
# cnts24only <- rowSums(mycounts24only$wide[, 1:2])
rm(mycounts)
# same as above for F4s
# cnts21only <- rowSums(mycounts21only$wide[, which(meta_sub$plant == "M82C")])
# cnts22only <- rowSums(mycounts22only$wide[, which(meta_sub$plant == "M82C")])
# cnts23only <- rowSums(mycounts23only$wide[, which(meta_sub$plant == "M82C")])
# cnts24only <- rowSums(mycounts24only$wide[, which(meta_sub$plant == "M82C")])
# nlibs <- nrow(meta_sub) # 88
# # save(mycounts, cnts21, cnts24, file = paste0(prefix, mygenome, "_counts.rdata"))
# df_counts <- data.frame(srnas21 = cnts21only,srnas22 = cnts22only,srnas23 = cnts23only,srnas24 = cnts24only)


# tmp <- cbind(slydf, df_counts)



# tmp2 <- tmp %>%
#   filter(!!sym(factor1) != "none") %>%
#   group_by(!!sym(factor1)) %>%
#   summarize(sum21 = sum(srnas21),
#             sum22 = sum(srnas22),
#             sum23 = sum(srnas23),
#             sum24 = sum(srnas24)
#             )

# }}}

#------------------------------------ total sRNA profile (mapped+unmapped) {{{
  #total
  srnas_tot <- slydf %>%
    summarise(across(contains("count"), ~ sum(.x))) %>%
    pivot_longer(
                 cols = everything(),
                 names_to = c("plant", "size", "type"),
                 names_pattern = "srnas_(.*)_(.*)_(.*)",
                 values_to = "count") %>%
    mutate(across(where(is.character), as.factor))

gg <- ggplot(srnas_tot, aes(color = plant, 
                       y = count,
                       group = plant,
                       x = size)) +
geom_point() +
geom_line() +
theme_bw() +
theme(legend.position = "right",
      axis.text.x = element_text(angle = 300, hjust = 0),
      axis.ticks = element_blank()) +
      # scale_fill_manual(values = brewer.pal(4,"PRGn"))  +
scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6, digits = 2),
                   breaks = scales::pretty_breaks(n = 4))
ggsave(paste0("srna_size_counts_profile_total_raw", ".pdf"),
       gg,
       height = 4,
       width = 4)

# tmp from srnas_processing.sh
gg <- ggplot(cnt_summarized, aes(color = plant, 
                       y = count,
                       group = plant,
                       x = size)) +
geom_point() +
geom_line() +
theme_bw() +
theme(legend.position = "right",
      axis.text.x = element_text(angle = 300, hjust = 0),
      axis.ticks = element_blank()) +
      # scale_fill_manual(values = brewer.pal(4,"PRGn"))  +
scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6, digits = 2),
                   breaks = scales::pretty_breaks(n = 4))
ggsave(paste0("srna_size_counts_profile_total_unmapped", ".pdf"), gg, height = 4, width = 4)

# }}}

#------------------------------------ looking at break downs for size classes {{{
# setwd(file.path(path_srnas, "analysis_bins/"))
myF4 ="P4042"
myF4 ="P1512"
infix = ""
# gtfilter = "00"; genotypename = "slyc"
# gtfilter = "11"; genotypename = "penn"
myprefix <- "srnas_"
suffix <- "_class2"
setwd(file.path(path_base, "comparisons_using_big_table/sRNA_genomewide/", genotypename, "/srna_sizes/"))


factor1  <- paste0(myprefix, myF4, infix, suffix)
factor1 <- "TEs_Order_renamed"
# mod for dcl2 class
factor1="srnas_dcl2_class2"
factor1="srnas_dcl2_tmv_class2"
factor1="srnas_M82_size_class"
# factor1="srnas_dcl2_tmv_25bp_count"

slydf$srnas_dcl2_WT_20bp_count <- slydf$srnas_dcl2_WT_20_class2
slydf$srnas_dcl2_WT_23bp_count <- slydf$srnas_dcl2_WT_23_class2
slydf$srnas_dcl2_WT_25bp_count <- slydf$srnas_dcl2_WT_25_class2
slydf$srnas_dcl2_d2Mutant_20bp_count <- slydf$srnas_dcl2_d2Mutant_20_class2
slydf$srnas_dcl2_d2Mutant_23bp_count <- slydf$srnas_dcl2_d2Mutant_23_class2
slydf$srnas_dcl2_d2Mutant_25bp_count <- slydf$srnas_dcl2_d2Mutant_25_class2
slydf$srnas_dcl2_WT_20_count <- NULL
slydf$srnas_dcl2_WT_23_count <- NULL
slydf$srnas_dcl2_WT_25_count <- NULL
slydf$srnas_dcl2_d2Mutant_20_count <- NULL
slydf$srnas_dcl2_d2Mutant_23_count <- NULL
slydf$srnas_dcl2_d2Mutant_25_count <- NULL

#------------------------------------ export srnas bins as table for supp data

all_F4s_tables <- list()
for (myF4 in F4s_srnas) {
  #   myF4 <- F4s_srnas[7]
  message(myF4)
  genotype <- sym(paste0("genotype_", myF4))
  factor1  <- paste0(myprefix, myF4, "", "_class")
  # factor1  <- paste0(myprefix, myF4, infix, "_class")
  #   rnaseq_a  <- paste0(mydata, myF4, "_class")
  #save number of DEG/DET for each F4
# colnames(mydf_filtera)
  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter) %>%
    filter(!!sym(paste0(factor1)) == "DESL")

  all_F4s_tables[[paste0("DESLs_M82_vs_", myF4)]] <- mydf_filtera %>%
  # tmp <- mydf_filtera %>%
    dplyr::select(c(1:28, c("srnas_dcl2_logFC", "srnas_dcl2_class2"), contains(myF4))) %>%
    dplyr::select(!(contains("meth") | ends_with("count") |contains("unique") | contains("rnaseq") | contains("class3") | ends_with("class") | contains("class_")))
}
map(all_F4s_tables, dim)

library(openxlsx)
openxlsx::write.xlsx(all_F4s_tables, "srnas_DESL_all_F4s.xlsx")




for (myF4 in F4s_srnas) {
  #   myF4 <- F4s_srnas[7]

  message(myF4)
  genotype <- sym(paste0("genotype_", myF4))
  factor1  <- paste0(myprefix, myF4, infix, suffix)
  # factor1  <- paste0(myprefix, myF4, infix, "_class")
  #   rnaseq_a  <- paste0(mydata, myF4, "_class")
  #save number of DEG/DET for each F4
  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter) %>%
    filter(!!sym(paste0(factor1)) != "inbetween")
    # filter(!!sym(paste0(factor2)) != "inbetween")

  message(nrow(mydf_filtera))

  # for (factor2 in factor2s) {

  # }
  #by factor1
  # table the counts sRNA species for each F4, category and size class 
   srnas_fac1 <- mydf_filtera %>%
  # srnas_fac1 <- slydf %>% # change for dcl2
    group_by(!!sym(factor1)) %>%
    summarise(across(contains("count"), ~ sum(.x)))
  # dplyr::select(1, contains("P3611"))
  # TEs_Order_renamed srnas_P3611_21bp… srnas_P3611_22bp… srnas_P3611_23b… srnas_P3611_24b… srnas_P3611_20b… srnas_P3611_25b…
  # <fct>                         <dbl>             <dbl>            <dbl>            <dbl>            <dbl>            <dbl>
  # 1 TIR                          186755            221728           393447           994157            99391            26162
  # 2 LTR                          968874           1530302          1029346          2311870           557211           167481

  srnas_fac1_long <- srnas_fac1 %>%
    pivot_longer(
                 cols = -1,
                 names_to = c("plant", "size", "type"),
                 names_pattern = "srnas_(.*)_(.*)_(.*)",
                 values_to = "count") %>%
  drop_na() %>%
  filter(!!sym(paste0(factor1)) != "inbetween") %>%
  # mod for no-non!!
  # filter(!!sym(paste0(factor1)) != "non-D2DESL") %>%
  # filter(!!sym(paste0(factor1)) != "non-DESL") %>%
  # filter(srnas_dcl2_class2 %in% c("-D2DESL", "+D2DESL")) %>%
   # filter(size %in% c("21bp", "22bp", "23bp", "24bp")) %>%
   # filter(plant %in% c("dcl2_WT", "dcl2_d2Mutant", "dcl2_tmv" )); libs="dcl2_libs" # change for dcl2
  # mod for dcl2
  # filter(plant %in% c("dcl2_WT", "dcl2_d2Mutant" )); libs="dcl2_libs" # change for dcl2
  filter(plant %in% c("M82C", myF4)); libs=paste0(myF4, "_libs")
  # ! if including Penn
# filter(plant %in% c("PenneC", "M82C", myF4)); libs=paste0(myF4, "_libs")
# TEs_Order_renamed plant size  type    count
# <fct>             <chr> <chr> <chr>   <dbl>
# 1 TIR               P3611 21bp  count  186755
# 2 TIR               P3611 22bp  count  221728


# plotting counts (TMM normalized)
gg <- ggplot(srnas_fac1_long, aes(y = count, 
                                  # ! if Penn included
                                  # color = fct_relevel(plant, "M82C", "PenneC")
                                  color = plant, 
                                  group = plant,
                                  x = size)) +
geom_point() +
geom_line() +
facet_grid(vars(), vars(!!sym(factor1))) +
gglayers2 +
# scale_y_continuous(trans='log10',
scale_y_continuous(
                   labels = scales::unit_format(unit = "k", scale = 1e-3, accuracy = 1)) +
                   # breaks = scales::pretty_breaks(n = 12)) +
labs(title = myF4, 
     x = "",
     y = "sRNA counts (normalized)")

ggsave(paste0("srna_size_counts_tmm_lineplot_by", factor1, "_", libs, ".pdf"), gg, height = 2.2, width = 5)
# ggsave(paste0("srna_size_counts_tmm_lineplot_by", factor1, "_", libs, ".pdf"), gg, height = 4, width = 7)
# ggsave(paste0("srna_size_counts_tmm_lineplot_by", factor1, "_", libs, "_no-non.pdf"), gg, height = 4, width = 7)
# ggsave(paste0("srna_size_counts_tmm_lineplot_by", factor1, "_", libs, "_non-non.pdf"), gg, height = 2.5, width = 4)
}

# barplot 
gg <- ggplot(srnas_fac1_long, aes(x = !!sym(factor1), 
                                  y = count,
                                  fill = size)) +
geom_bar(stat="identity", position="fill") +
facet_grid(. ~ plant) +
gglayers2 +
scale_y_continuous(labels = scales::percent_format(scale =100, acc = 1)) +
labs(title = myF4,
     x = "",
     y = "sRNA size proportion")
# ggsave(paste0("srna_size_counts_barplot_by", factor1,"_",libs, ".pdf"), gg, height = 4, width = 4)


# # normalizing within plant and factor1 (e.g. DESL)
# sfl_norm <- srnas_fac1_long %>%
#   group_by(!!sym(factor1), plant) %>%
#   # mutate(prop = scale(count))
#   mutate(!!sym(paste0(factor1, "_freq")) := round(100 * count / sum(count),2))

# gg <- ggplot(sfl_norm, aes(y = !!sym(paste0(factor1, "_freq")), 
#                            color = plant,
#                            group = plant,
#                            x = size)) +
# geom_point() +
# geom_line() +
# # facet_grid(. ~ !!sym(factor1)) +
# # facet_grid(. ~ !!sym(factor1)) +
# facet_grid(vars(), vars(!!sym(factor1))) +
# # facet_grid(. ~ srnas_P1512_class2) +
# gglayers3 +
# labs(title = myF4, 
#      x = "",
#      y = "sRNA size proportion")
# ggsave(paste0("srna_size_counts_lineplot_by", factor1, ".pdf"), gg, height = 3, width = 5)

# }}}

#------------------------------------ 2d print barplot and percentage plots {{{
# 2D normalizing within plant and factor1 (e.g. DESL) and annotation
# factor1 <- "TEs_Order_renamed"
factor1="srnas_dcl2_class2"
factor2s <- c("TEs_Order_renamed", "srnas_dcl2_class2", "TEs_Superfamily", "rnaseq_dcl2_class2", "srnas_M82_size_class", "annotation", "chromatin_state", "alltrue")
factor2 <- c("TEs_Superfamily")
factor2 <- c("TEs_Order_renamed")
factor2 <- c("srnas_dcl2_class2")
factor2 <- c("srnas_tmv_class2")
factor2s <- factor2s[factor1 != factor2s]

mystats <- NULL
mystatsbin <- NULL
for (myF4 in F4s_srnas) {
  #   myF4 <- F4s_srnas[7]

  message(myF4)
  genotype <- sym(paste0("genotype_", myF4))
  factor1  <- paste0(myprefix, myF4, infix, suffix)
  # factor1  <- paste0(myprefix, myF4, infix, "_class")
  #   rnaseq_a  <- paste0(mydata, myF4, "_class")
  #save number of DEG/DET for each F4
  mydf_filtera <- slydf %>%
    filter(!!genotype == gtfilter) %>%
    filter(!!sym(paste0(factor1)) != "inbetween") %>%
    filter(!!sym(paste0(factor2)) != "inbetween")


  tmp <- mydf_filtera %>%
    filter(!!sym(factor2) == "-D2DESL") %>%
    dplyr::count(!!sym(paste0(factor1)), TEs_Order_renamed) %>%
    group_by(!!sym(paste0(factor1))) %>%
    mutate(n_freq = n/sum(n)) %>%
    mutate(plant = myF4) %>%
    rename(F4_class = !!sym(paste0(factor1)))
   mystatsbin <- rbind(mystatsbin, tmp)
# }

  message(nrow(mydf_filtera))

  # for (factor2 in factor2s) {

    # with(mydf_filtera, table_prop(TEs_Order_renamed, srnas_P4042_class2, mymargin = 2))
    with(mydf_filtera, table(TEs_Order_renamed, srnas_P4042_class2))
    with(mydf_filtera, table(TEs_Order_renamed, srnas_P4042_class2, srnas_M82_size_class))
    with(mydf_filtera, table(TEs_Superfamily, srnas_P4042_class2))

    srnas_fac1 <- mydf_filtera %>%
      # srnas_fac1 <- slydf %>% # for dcl2
      group_by(!!sym(factor1), !!sym(factor2)) %>%
      summarise(across(contains("count"), ~ sum(.x))) %>%
      ungroup()


    # for fig3 plot (panel D)

    srnas_fac1_eprv <- mydf_filtera %>%
      # srnas_fac1 <- slydf %>% # for dcl2
      filter(TEs_Order_renamed == "Pararetrovirus")  %>%
      group_by(!!sym(factor1), !!sym(factor2)) %>%
      summarise(across(contains("count"), ~ sum(.x))) %>%
      ungroup()

    srnas_fac1_filt <- srnas_fac1 %>%
      filter(!!sym(factor2) == "-D2DESL") %>%
      select(!!sym(paste0(myprefix, myF4, infix, "_22bp_count")), !!sym(paste0(myprefix, "M82C", infix, "_22bp_count")))
    srnas_fac1_eprv_filt <- srnas_fac1_eprv %>%
      filter(!!sym(factor2) == "-D2DESL") %>%
      select(!!sym(paste0(myprefix, myF4, infix, "_22bp_count")), !!sym(paste0(myprefix, "M82C", infix, "_22bp_count")))
    tmp <- srnas_fac1_eprv_filt/srnas_fac1_filt
    # tmp[,3] <- srnas_fac1_eprv_filt[,1]
    colnames(tmp) <- c("F4", "M82")
    tmp[,"class2"] <- c("-DESL", "+DESL", "non-DESL")
    tmp[,"plant"] <- myF4

    mystats <- rbind(mystats, tmp)

    srnas_fac1_long_eprv <- srnas_fac1_eprv %>%
      pivot_longer(
                   cols = -c(1:2),
                   names_to = c("plant", "size", "type"),
                   names_pattern = "srnas_(.*)_(.*)_(.*)",
                   values_to = "count") %>%
    drop_na() %>%
    # filter(plant %in% c("PenneC", "M82C", myF4))
    filter(plant %in% c("M82C", myF4)) %>%
    # mutate(plant = paste0(plant, "_eprv")) %>%
    mutate(myfilter = "EPRV")


  srnas_fac1_long <- srnas_fac1 %>%
    pivot_longer(
                 cols = -c(1:2),
                 names_to = c("plant", "size", "type"),
                 names_pattern = "srnas_(.*)_(.*)_(.*)",
                 values_to = "count") %>%
  drop_na() %>%
  mutate(myfilter = "All") %>%
  # filter(size %in% c("21bp", "22bp", "23bp", "24bp")) %>%
  # filter(plant %in% c("PenneC", "M82C", myF4))
   filter(plant %in% c("M82C", myF4)); libs=paste0(myF4, "_libs")
  # filter(plant %in% c("dcl2_WT", "dcl2_d2Mutant", "dcl2_tmv" )); libs="dcl2_libs" # change for dcl2
# filter(plant %in% c("dcl2_WT", "dcl2_d2Mutant" )); libs="dcl2_libs"
# srnas_fac1_long$plant
# srnas_fac1_long <- rbind(srnas_fac1_long, srnas_fac1_long_eprv)

# gg <- ggplot(srnas_fac1_long, aes(x = !!sym(factor1), 
#                                   y = count,
#                                   fill = size)) +
# geom_bar(stat="identity", position="fill") +
# facet_grid(vars(!!sym(factor2)), vars(plant)) +
# gglayers2 +
# scale_y_continuous(labels = scales::percent_format(scale =100, acc = 1)) +
# labs(title = myF4,
#      x = "",
#      y = "sRNA size proportion")
# ggsave(paste0("srna_size_counts_barplot_by_", factor2, "_", factor1, libs, ".pdf"), gg, height = 10, width = 4)

# tmm normalized counts
gg <- ggplot(srnas_fac1_long, aes(y = count, 
                                  color = plant,
                                  group = plant,
                                  x = size)) +
geom_vline(xintercept = c("21bp", "22bp", "24bp"), color = "grey") + # vertical line +
geom_point(size = 1) +
# geom_line(aes(linetype=myfilter)) +
geom_line() +
geom_line(data = srnas_fac1_long_eprv, linetype = "dotted") +
# geom_line(data = srnas_fac1_long_eprv, aes(linetype = "dotted")) +
# facet_grid(. ~ !!sym(factor1)) +
# facet_grid(. ~ !!sym(factor1)) +
facet_grid(vars(!!sym(factor2)), vars(!!sym(factor1)), scale = "free") +
# facet_grid(vars(!!sym(factor1)), vars(!!sym(factor2)), scale = "free") +
gglayers2 +
# scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3, digits = 0),
#                    breaks = scales::pretty_breaks(n = 4)) +
# scale_y_continuous(trans='log10',
scale_y_continuous(
                   labels = scales::unit_format(unit = "k", scale = 1e-3, accuracy = 1)) +
labs(title = myF4,
     x = "",
     y = "sRNA counts (normalized)")
# ggsave(paste0("srna_size_counts_tmm_lineplot_by_", factor2, "_", factor1, libs, ".pdf"), gg, height = 8, width = 6)
# ggsave(paste0("srna_size_counts_tmm_lineplot_by_", factor2, "_", factor1, libs, "_log.pdf"), gg, height = 8, width = 6)

ggsave(paste0("srna_size_counts_tmm_lineplot_by_", factor2, "_", factor1, libs, ".pdf"), gg, height = 4.5, width = 5)
# ggsave(paste0("srna_size_counts_tmm_lineplot_by_", factor2, "_", factor1, ".pdf"), gg, height = 5, width = 15)

  }

write_csv(mystats, "fraction_eprv_total_22bp_-D2SL.csv")
write_csv(mystatsbin, "fraction_eprv_total_bins_-D2SL.csv")

mystats_long <- mystats %>%
  pivot_longer(cols = 1:2, names_to = "condition", values_to = "fraction")

gg <- ggplot(mystats_long, aes(x = class2, y = fraction*100, shape = plant, group = class2)) +
  gglayers2 +
                scale_shape_manual(values = F4_shapes) +
                scale_fill_manual(values = TRIBE_colors) +
  # facet_grid(. ~ condition) +
  # stat_summary(geom = "bar", fill = "lightgrey") +
  geom_boxplot(aes(fill = class2,alpha = 0.7)) +
  # stat_summary(geom = "bar") +
  geom_point() +
  labs(title = "",
     x = "",
     y = "% EPRV of 22nt sRNAs for -D2SL bins")
ggsave(paste0("fraction_eprv_total_22bp_-D2SL.pdf"), gg, height = 4, width = 3)

gg <- ggplot(filter(mystatsbin, TEs_Order_renamed == "Pararetrovirus"),
             aes(x = F4_class, y = 100*n_freq, shape = plant, group = F4_class)) +
  gglayers2 +
                scale_shape_manual(values = F4_shapes) +
                scale_fill_manual(values = TRIBE_colors) +
  # facet_grid(. ~ condition) +
                geom_boxplot(aes(fill = F4_class,alpha = 0.7)) +
  # stat_summary(geom = "bar", fill = "lightgrey") +
  # stat_summary(geom = "bar") +
  geom_point() +
  labs(title = "",
     x = "",
     y = "% EPRV bins for -D2SL")
ggsave(paste0("fraction_eprv_total_bins_-D2SL.pdf"), gg, height = 4, width = 3)

# counting eprv proportion using bin counts


   dplyr::count(srnas_P4042_class2, TEs_Order )

# sfl_norm <- srnas_fac1_long %>%
#   # group_by(!!sym(factor1), plant, dcl2_class) %>%
#   group_by(!!sym(factor1), plant, !!sym(factor2)) %>%
#   mutate(!!sym(paste0(factor1, "_freq")) := round(100 * count / sum(count),2))

# gg <- ggplot(sfl_norm, aes(y = !!sym(paste0(factor1, "_freq")), 
#                            color = plant,
#                            group = plant,
#                            x = size)) +
# geom_point() +
# geom_line() +
# # facet_grid(. ~ !!sym(factor1)) +
# # facet_grid(. ~ !!sym(factor1)) +
# facet_grid(vars(!!sym(factor2)), vars(!!sym(factor1))) +
# # facet_grid(. ~ srnas_P1512_class2) +
# theme_bw() +
# theme(legend.position = "right",
#       axis.text.x = element_text(angle = 300, hjust = 0),
#       axis.ticks = element_blank()) +
# scale_y_continuous(labels = scales::percent_format(scale =1, acc = 1)) +
# labs(title = myF4,
#      x = "",
#      y = "sRNA size proportion")
# ggsave(paste0("srna_size_counts_lineplot_by_", factor2, "_", factor1, ".pdf"), gg, height = 8, width = 5)


# }}}

#------------------------------------ old way: print barplot and perc plot for sizes {{{

tmp2 <- srnas_fac1
tmp3 <- round(prop.table(as.matrix(tmp2[,2:5]), margin = 1)*100,1)
rownames(tmp3) <- tmp2 %>% pull(1)
write_csv(as.data.frame(tmp3), path = paste0("srna_size_counts_", factor1, ".csv"))
 print(tmp3)
               # sum21 sum22 sum23 sum24
# DNA              7.5  10.2  20.4  61.9
# LTR             11.2  16.3  18.4  54.1
# none            16.9  12.7  17.2  53.2
# nonLTR           7.3   7.9  19.8  65.0
# Pararetrovirus  11.7  64.3   6.8  17.3
# PHG              7.9   7.1  19.7  65.3
# Rolling         75.4  10.6   3.7  10.3
# SSR              6.5   7.3  21.8  64.4
# TE               7.8   7.0  20.1  65.1
# Unknown          6.2   9.6  21.0  63.2
# <NA>            13.2  18.5  17.4  50.9
# <NA>            12.2  13.6  18.3  55.8

library(tidyr)
tmp4 <- as.data.frame(tmp3) %>%
  tibble::rownames_to_column() %>%
  gather(size_class, proportion, -rowname) %>%
  mutate(size_class = str_replace(size_class, "sum",""))
  # mutate(TEs_Order_renamed = fct_relevel(TEs_Order_renamed, c("Pararetrovirus", "Helitron", "Other_TEs"), after = Inf))

gg <- ggplot(tmp4, aes(fill = size_class, 
                       y = proportion,
                       x = rowname)) +
# geom_point() +
geom_bar(stat="identity", position="fill") + 
# geom_line() +
theme_bw()+
theme(legend.position = "right",
      axis.text.x = element_text(angle = 300, hjust = 0),
      axis.ticks = element_blank()) +
      scale_fill_manual(values = brewer.pal(4,"PRGn"))  +
scale_y_continuous(labels = scales::percent,
                   breaks = scales::pretty_breaks(n = 5))

# ggsave("srna_size_counts_TEs_order.pdf",
ggsave(paste0("srna_size_counts_", factor1, ".pdf"),
       gg,
       height = 4,
       width = 4)

# testing 21/24 ratio
# setwd(path_map)
# tmp <- readGAlignments(fls[1])

# bin.counts.wide <- matrix(NA, nr = length(mybins), 
#                               nc = nlibs)
# colnames(cnts[,1:nlibs]) <- meta$lib
# rownames(bin.counts.wide) <- mybins$ID
# # filter counts

cnts <- cnts %>%
  mutate(myfilter = rowSums(dplyr::select(., contains("bam"))) > 50)
sly_bins$myfilter <- cnts$myfilter

tmp <- meta %>% 
  filter(generation == "F4")

# }}}

#------------------------------------ read in idividual srnas form edgeR {{{
#dcl2
load(file="/projects/TRIBE/srnas/dcl2_zhengming/cnt_dcl2.rdata") #cnt_dcl2

lib_sizes <- colSums(cnt_dcl2[,3:6])

cnt_dcl2_filt <- cnt_dcl2 %>%
  filter(width %in% c(21,22)) %>%
  filter(rowSums(dplyr::select(., contains("-"))) > 100)

  group <- c("dcl","dcl","wt","wt")

  dgList_tmp <- DGEList(counts=cnt_dcl2[,3:6],
                    group=group,
                    remove.zeros = F,
                    genes=cnt_dcl2$sequence
  )
  dgList_tmp <- calcNormFactors(dgList_tmp, method="TMM")
  mysamples <- dgList$samples

  dgList <- DGEList(counts=cnt_dcl2_filt[,3:6],
                    group=group,
                    remove.zeros = F,
                    genes=cnt_dcl2_filt$sequence
  )
  dgList$samples <- mysamples

prefix <- "sRNA_species_filter_"
design <- model.matrix(~group)
dgList <- estimateDisp(dgList, design)
png(file=paste0(prefix, "dcl2", "_BCV", ".png"), width=1000, height=1000, type = "cairo-png")
plotBCV(dgList)
dev.off()
fit <- glmFit(dgList,design)
lrt_dcl2 <- glmLRT(fit, coef=2)
deGenes <- decideTestsDGE(lrt_dcl2, p=0.00001)
deGenes <- rownames(lrt_dcl2)[as.logical(deGenes)]
png(file=paste0(prefix, "dcl2", "DE", ".png"), width=1000, height=1000, type = "cairo-png")
plotSmear(lrt_dcl2, de.tags=deGenes)
dev.off()

srnas_species_dcl2 <- lrt_dcl2$table %>%
    mutate(FDR = p.adjust(PValue, method = "BH")) %>%
    mutate(updown = ifelse(logFC < 0, "-", "+")) %>%
    mutate(logCPM = round(logCPM,2)) %>%
    mutate(logFC = round(logFC,2))

srnas_species_dcl2_filt <- cbind(srnas_species_dcl2, cnt_dcl2_filt) %>%
  filter(FDR < 1e-8) %>%
  filter(updown == "+")

tmp3 <- DNAStringSet(table_all$sequence)
names(tmp3) <- paste0("dcl2_-", DNAStringSet(table_all$sequence))
writeXStringSet(tmp3, "dcl2-_srnas_species_stringent.fa", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

#}}}

#------------------------------------ read inP4042 individual srnas {{{
load(file.path(path_srnas, "cnt_tribe.rdata"), envir = myenv <- new.env()) #cnt
ls(myenv)
str(myenv$cnt)

lib_sizes <- myenv$cnt %>%
  summarise(across(where(is.numeric), ~ sum(.x))) %>%
  unlist()

cnt_tribe_filt <- cnt %>%
  filter(width %in% c(21,22)) %>%
  filter(rowSums(dplyr::select(., contains("-"))) > 100)

#------------------------------------ annotatin sRNA species with TEs etc
library(GenomicAlignments)
p1 <- ScanBamParam( what=scanBamWhat())
myfile <- paste0("all", "_srnas_species_100.fa")
myfile2 <- file.path(path_srnas, "mirnas_species", myfile)
# [1] "/projects/TRIBE/srnas/P4042-_srnas_species_stringent_500.fa"

tmp3 <- DNAStringSet(cnt_tribe_filt$sequence)
names(tmp3) <- paste0(DNAStringSet(cnt_tribe_filt$sequence))

writeXStringSet(tmp3, file = myfile2, append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

lyc="/data/public_data/tomato/S_lycopersicum/S_lycopersicum_chromosomes.3.00"
system( paste0("bowtie --wrapper basic-0 -v 0 -p 24 -k 20 -f -S ", 
              lyc, " ",
              myfile, " ", myfile,".sam" ))
# reads with at least one reported alignment: 51721 (63.06%)
# reads that failed to align: 30293 (36.94%)
# Reported 401050 alignments
system( paste0("samtools view -b ", myfile, ".sam > ", myfile, "_multi.bam"))
readin <- readGAlignments(paste0(myfile,"_multi.bam"), param=p1)

readinm <- readin[strand(readin) == "-",]
readinp <- readin[strand(readin) == "+",]
readinm@elementMetadata$seq <- reverseComplement(readinm@elementMetadata$seq)
readin <- c(readinm,readinp)

readin <- overlap_annot(readin, sly_bins, "TEs_Order", "TEs_Order")
readin <- overlap_annot(readin, sly_bins, paste0("srnas_",myF4,"_class2"), paste0("srnas_",myF4,"_class2"))
readin <- overlap_annot(readin, sly_bins, "srnas_dcl2_class2", "srnas_dcl2_class2")
readin <- overlap_annot(readin, sly_bins, "genotype_P4042", "genotype_P4042")

rdf <- as.data.frame(readin) %>%
  distinct(seq, TEs_Order, .keep_all= TRUE) %>%
  dplyr::rename(sequence = seq)

tmp <- cnt_tribe_filt %>%
  left_join(rdf, by = "sequence")
  # filter(srnas_P4042_class2 == "+DESL")
  # dplyr::count(srnas_P4042_class2, TEs_Order )

# adding DE for M82/Penn vs P4042
counts <- tmp %>%
  dplyr::select(contains("Penn") | contains("M82") | contains("P4042-")) %>%
  mutate(across(everything(), as.numeric))

groups <- factor(c(rep("Parent", 4), rep("P4042", 2)))
y <- DGEList(counts=counts,group=groups)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
colnames(lrt$table) <- paste("Parent_vs_P4042", colnames(lrt$table), sep="_")

tmp2 <- cbind(tmp, lrt$table)

save(cnt_tribe_filt, file=file.path(path_srnas, "counts_raw_srnas_tribe.rdata")) #cnt
# load(file=file.path(path_srnas, "counts_raw_srnas_tribe.rdata")) #cnt
write_csv(tmp2, path = paste0("all_srnas_", mygenome, "_", "_srnas_species_", "all", "_filt_table_annotation.csv"))

#segmentSeq
# mirRNAs_path <- "/data/public_data/tomato/additional_resources/miRNAs_lyc_and_penn/miRNAs_sly_spenn.fa"
# mirRNAs_path <- "/data/public_data/mirBase/mature_miRNAs_r22_UtoT.fa"
mirRNAs_path <- "/data/public_data/mirBase/mature_miRNAs_r22_UtoT_Sly.fa"
# sly-miR6026
mirnas <- readDNAStringSet(mirRNAs_path)

#only flagging
cnt_tribe_filt$ismiR <- as.character(cnt_tribe_filt$sequence) %in% as.character(mirnas)
cnt_tribe_filt$miRNA <- rep(NA, nrow(cnt_tribe_filt)) 

idx <- which(cnt_tribe_filt$ismiR)

#annotating sRNAs with miRNA name if found
#idea: optimizing by which statement of found vector and back ref..
for (i in idx) {
	cnt_tribe_filt$miRNA[i] <- names(as.character(mirnas)[as.character(mirnas) %in% cnt_tribe_filt$sequence[i]])[1]
}

#------------------------------------ sRNA profile for all F4s

  idx2 <- grep("-", colnames(cnt_tribe_filt))
  cnt_tribe_myF4_filte <- cnt_tribe_filt[,c(idx2)]
  # cnt_tribe_myF4_filte <- cnt_tribe_filt[,c(1:3, idx)]
  group <- c("Parent","Parent","Parent","Parent", rep("F4",18))

  dgList_tmp <- DGEList(counts=cnt_tribe_myF4_filte,
                    group=group,
                    remove.zeros = F,
                    genes=cnt_tribe_filt$sequence
  )

  # replacing pre-calculated lib sizes
  tmp <- dgList_tmp$samples
  tmp$lib.size <- lib_sizes[rownames(tmp)]
  dgList_tmp$samples <- tmp

meta_sub %>%
  filter(generation == "F4")

cntsnc  =  cpm(dgList_tmp, normalized.lib.sizes = TRUE)
group <- str_replace(colnames(cntsnc),"-.*", "")
group <- str_replace(group,".*_", "")

srnas_species <- cbind(as.data.frame(round(cntsnc,3)), miRNA = cnt_tribe_filt$miRNA)

tmp <- srnas_species %>%
  filter(!is.na(miRNA))

write_csv(tmp, path = paste0("mirnas", mygenome, "_", "_srnas_species_", "all", "_filt_table.csv"))
read_csv(paste0("mirnas", mygenome, "_", "_srnas_species_", "all", "_filt_table.csv"))

pdf(file="mirnas.pdf", width = 8, height = 5) # or other device
plot(unlist(tmp[1,1:22]))
tmp2=tmp[1,]
ggplot(tmp[1,1:22], aes(y=))
plot((tmp[1,1:22]))
dev.off()

#------------------------------------ subsetting F4s

  myF4 <- "P4042"
  myF4 <- "P4041"
  myF4 <- "P3611"
  myF4 <- "P2562"
  myF4 <- "P2561"
  myF4 <- "P1512"
  myF4 <- "P3612"

for (myF4 in F4s_srnas) {
  group <- c("Parent","Parent","Parent","Parent", myF4, myF4)
  # group <- c("Parent","Parent", myF4, myF4)

  idx2 <- grep(myF4, colnames(cnt_tribe_filt))
  cnt_tribe_myF4_filte <- cnt_tribe_filt[,c(1:5, idx2)]
  # cnt_tribe_myF4_filte <- cnt_tribe_filt[,c(1:3, idx)]

  dgList_tmp <- DGEList(counts=cnt_tribe_myF4_filte[,2:7],
  # dgList_tmp <- DGEList(counts=cnt_tribe_myF4_filte[,2:5],
                    group=group,
                    remove.zeros = F,
                    genes=cnt_tribe_filt$sequence
  )
  # replacing pre-calculated lib sizes
  tmp <- dgList_tmp$samples
  tmp$lib.size <- lib_sizes[rownames(tmp)]
  dgList_tmp$samples <- tmp


prefix <- paste0("sRNA_species_filter_libcorrected_", myF4)
dgList <- calcNormFactors(dgList_tmp, method="TMM")
design <- model.matrix(~group)
dgList <- estimateDisp(dgList, design)


png(file=paste0(prefix, "_BCV", ".png"), width=1000, height=1000, type = "cairo-png")
# plotBCV(dgList)
plotBCV(dgList[idx,])
dev.off()
fit <- glmFit(dgList,design)
lrt_myF4 <- glmLRT(fit, coef=2)

srnas_species_myF4 <- cbind(lrt_myF4$table, cnt_tribe_myF4_filte, miRNA = cnt_tribe_filt$miRNA) %>%
    mutate(logFC = -round(logFC,2)) %>%
    mutate(FDR = p.adjust(PValue, method = "BH")) %>%
    mutate(updown = ifelse(logFC < 0, "-", "+")) %>%
    mutate(logCPM = round(logCPM,2))

srnas_species_myF4_filt <- srnas_species_myF4 %>%
  filter(FDR < 5e-1) %>%
  filter(updown == "+") %>%
  filter(rowSums(dplyr::select(., contains("M82C-"))) > 300) # filter only high M82!
  # logFC logCPM       LR      PValue       FDR updown              sequence 01_KM82C-1 01_KM82C-2 02_KPenneC-1 02_KPenneC-2 22_KP4042-1 22_KP4042-2
# 1 10.40   5.37 5.430143 0.019792056 0.3038507      + AAAAAGCTATAGGTCAATTGA        252        292            0            0           0           0

srnas_species_myF4_filt_mir<- srnas_species_myF4 %>% 
  filter(FDR < 5e-1) %>%
  filter(!is.na(miRNA))

write_csv(srnas_species_myF4_filt, path = paste0(prefix, mygenome, "_", "_srnas_species_", myF4, "_filt_table.csv"))
write_csv(srnas_species_myF4_filt_mir, path = paste0(prefix, mygenome, "_", "_srnas_species_", myF4, "_filt_table_mir.csv"))
# [1] "sRNA_species_filter_P4042SolLyc_srnas_species_P4042_filt_table.csv"
save(srnas_species_myF4_filt, file = file.path(path_srnas, paste0(prefix, mygenome, "_srnas_species_",myF4, "_filt")))
# load("/projects/TRIBE/srnas/sRNA_species_filter_P4042SolLyc_srnas_species_P4042_filt")
}

myfile <- paste0(myF4, "-_srnas_species_stringent_300.fa")
myfile2 <- file.path(path_srnas, myfile)
# [1] "/projects/TRIBE/srnas/P4042-_srnas_species_stringent_500.fa"

tmp3 <- DNAStringSet(srnas_species_myF4_filt$sequence)
names(tmp3) <- paste0(myF4, "_-", DNAStringSet(srnas_species_myF4_filt$sequence))

writeXStringSet(tmp3, file = myfile2, append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

  # tags <- topTags(lrt_dcl2, n = 1e5, p.value = 1e-9)
  # write_csv(as.data.frame(tags), path = paste0(prefix, mygenome, "_", myF4, "_table.csv"))

# load(file = paste0(prefix, mygenome, "_bins_raw_processed.rdata"))
# load(file = file.path(path_srnas, "analysis_bins" , paste0(prefix, mygenome, "_bins_raw_processed_unique.rdata")))

# }}}
save(srna_cnts_bin, file = paste0(prefix, mygenome, "_bins_raw_processed.rdata"))
save(lrt_dcl2, mycounts_dcl2, file = paste0(prefix, mygenome, "_bins_raw_processed_dcl2.rdata"))
#}}}

# Are 21 or 22
# 	Are part of a -DESL in P4042
# 	Are part of a -D2SL
# 	The actual sRNA is going down in P4042
# 	The actual sRNA is going down in dcl2


#------------------- target analysis srnas {{{
#------------------------------------ load fa file and correlate it back with slydf

# merged=/data/public_data/tomato/merged_genomes/genomes_SL30_penn_organelles_merged
lyc="/data/public_data/tomato/S_lycopersicum/S_lycopersicum_chromosomes.3.00"

setwd(path_srnas)
system( paste0("bowtie --wrapper basic-0 -v 0 -p 24 -k 20 -f -S ", 
              lyc, " ",
              myfile, " ", myfile,".sam" ))
system( paste0("samtools view -b ", myfile, ".sam > ", myfile, "_multi.bam"))

system( paste0("bowtie --wrapper basic-0 -v 0 -m 1 -p 24 -k 1 -f -S ", 
              lyc, " ",
              myfile, " ", myfile,"_unique.sam" ))
system( paste0("samtools view -b ", myfile, "_unique.sam > ", myfile, "_unique.bam"))

library(GenomicAlignments)
p1 <- ScanBamParam( what=scanBamWhat())
# p1 <- ScanBamParam( what="seq")
readin <- readGAlignments(paste0(myfile,"_multi.bam"), param=p1)
readin <- readGAlignments(paste0(myfile,"_unique.bam"), param=p1)

readinm <- readin[strand(readin) == "-",]
readinp <- readin[strand(readin) == "+",]
readinm@elementMetadata$seq <- reverseComplement(readinm@elementMetadata$seq)
readin <- c(readinm,readinp)

# load(file = file.path(path_comparison, paste0("SolLyc", "_bins_plus_statsv23.rdata")))


sly_bins$srnas_dcl2_class2  <- slydf$srnas_dcl2_class2
sly_bins$srnas_P3611_class2 <- slydf$srnas_P3611_class2
sly_bins$srnas_P3612_class2 <- slydf$srnas_P3612_class2
sly_bins$srnas_P2562_class2 <- slydf$srnas_P2562_class2
sly_bins$srnas_P4041_class2 <- slydf$srnas_P4041_class2

sly_bins$srnas_P2561_class2 <- slydf$srnas_P2561_class2

readin <- overlap_annot(readin, sly_bins, "TEs_Order", "TEs_Order")
readin <- overlap_annot(readin, sly_bins, paste0("srnas_",myF4,"_class2"), paste0("srnas_",myF4,"_class2"))
readin <- overlap_annot(readin, sly_bins, "srnas_dcl2_class2", "srnas_dcl2_class2")

readin2 <- readin %>%
  as.data.frame %>%
  # dplyr::count(srnas_P4042_class2)
  dplyr::select(-c(strand.1, qwidth.1, qname, pos, mrnm, mpos, isize, qual, cigar, njunc, flag, rname)) %>%
  # filter(srnas_P2562_class2 == "-DESL") %>%
  filter(!!sym(paste0("srnas_",myF4,"_class2")) == "-DESL") %>%
  filter(srnas_dcl2_class2 == "-D2DESL") %>%
  distinct(seq, .keep_all= TRUE)
  
    # seqnames strand qwidth    start      end width strand.1      pos qwidth.1 mapq cigar.1 mrnm mpos isize                    seq      TEs_Order srnas_P4042_class2 srnas_dcl2_class2
# 1  SL3.0ch07      +     22 44644140 44644161    22        + 44644140       22  255     22M <NA>   NA     0 AAGTGTCTGATTTATCAAACTC           none              -DESL           -D2DESL
# 2  SL3.0ch01      +     22 76990255 76990276    22        + 76990255       22  255     22M <NA>   NA     0 ATGACGATCTTCTGGTCTGTGA Pararetrovirus              -DESL           -D2DESL
# 3  SL3.0ch06      -     22 49177916 49177937    22        - 49177916       22  255     22M <NA>   NA     0 CCACTTCTGGAGTCCACTTATG Pararetrovirus              -DESL           -D2DESL
# 4  SL3.0ch06      -     22  2280450  2280471    22        -  2280450       22  255     22M <NA>   NA     0 TATCGCAGAGTTTGATAAATCA Pararetrovirus              -DESL           -D2DESL
# 5  SL3.0ch09      +     22 68593406 68593427    22        + 68593406       22  255     22M <NA>   NA     0 TGGTCTAGAATCTTGAATTTCT            PHG              -DESL           -D2DESL
# 6  SL3.0ch05      +     22 37924953 37924974    22        + 37924953       22  255     22M <NA>   NA     0 TTAATGATCTCAATTGTAAATG Pararetrovirus              -DESL           -D2DESL
# 7  SL3.0ch06      -     22  2280662  2280683    22        -  2280662       22  255     22M <NA>   NA     0 TAAACAGTATTCAAGATTCGAA Pararetrovirus              -DESL           -D2DESL
# 8  SL3.0ch08      +     22 11885232 11885253    22        + 11885232       22  255     22M <NA>   NA     0 TTTAGGATCTAGAATAAAATTT            LTR              -DESL           -D2DESL
# 9  SL3.0ch01      -     22 76989577 76989598    22        - 76989577       22  255     22M <NA>   NA     0 GAAGAGCCCACAAGTTTATAAA Pararetrovirus              -DESL           -D2DESL
# 10 SL3.0ch11      -     21  2721390  2721410    21        -  2721390       21  255     21M <NA>   NA     0  TGTTCTACTCAAGATCGCAAA           none              -DESL           -D2DESL


table(readin2$srnas_dcl2_class2)

tmp3 <- DNAStringSet(readin2$seq)
names(tmp3) <- paste0("srnas_species_",myF4,"-DESL_dcl2_-", (readin2$seq), "_", (readin2$TEs_Order))
writeXStringSet(tmp3, paste0("srnas_species_",myF4,"-DESL_dcl2_verystringent_multi.fa"), append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
# go into Cleavland analysis with fa file
#--- target analysis srnas }}}


#------------------------------------ read in multiple/unique srnas form edgeR {{{
srna_cnts_bin <- srna_cnts_unique_bin <- list()

meta_sub <- meta %>%
  filter(include == "yes")
# filter(plant == "P3611") %>%
meta_sub2 <- meta_sub %>%
  dplyr::mutate(sub = mapped == "SolLyc"  )


cnts <- counts$wide
# read in all at once to determine libScalefactor

  # DE analysis for sRNA bins
  group <- subset(meta_sub2, sub)$plant

  # meta_sub$sub <- TRUE
  group <- meta_sub2$plant
  dgList <- DGEList(counts=cnts[,which(meta_sub2$sub)],
                    group=group[which(meta_sub2$sub)],
                    remove.zeros = F,
                    genes=cnts$ID
  )

  dgList <- calcNormFactors(dgList, method="TMM")
  lib.sizes <- dgList$sample$lib.size * dgList$samples$norm.factors
  # names(lib.sizes) <- rownames(dgList$sample)
  meta_sub$lib.sizes <- lib.sizes

cntsnc  =  cpm(dgList, normalized.lib.sizes = TRUE)
cnt_filt <- cntsnc[cnts$genome=="SL3.0-Genome",]

cnt_filt <- as.data.frame((cnt_filt)) %>%
  rename_with(~str_replace(.x,"_S.*", "_cpm")) %>%
  rename_with(~str_replace(.x,".*_K","srnas_"))

save(cnt_filt, file=file.path(path_srnas, "counts_cpm_alllineages_srnas_tribe.rdata")) # cnt_filt




myparent <- "PenneC"
# myparent <- "M82C"
# multiple
for (myF4 in F4s_srnas) {
  message(myF4)
  # myF4 <-  myF4s[7]
  # select only parents and F4s for time being:
  meta_sub2 <- meta_sub %>%
    dplyr::mutate(sub = mapped == "SolLyc" & (plant == myparent | plant == myF4))


  # DE analysis for sRNA bins
  group <- subset(meta_sub2, sub)$plant
  dgList <- DGEList(counts=cnts[,which(meta_sub2$sub)],
                    group=group,
                    remove.zeros = F,
                    genes=cnts$ID
  )

  dgList <- calcNormFactors(dgList, method="TMM")
  design <- model.matrix(~group)
  dgList <- estimateDisp(dgList, design)
  png(file=paste0(prefix, myF4, myparent, "_BCV", ".png"), width=1000, height=1000, type = "cairo-png")
  plotBCV(dgList)
  dev.off()
  fit <- glmFit(dgList,design)
  lrt <- glmLRT(fit, coef=2)
  deGenes <- decideTestsDGE(lrt, p=0.001)
  deGenes <- rownames(lrt)[as.logical(deGenes)]
  png(file=paste0(prefix, myF4, myparent, "DE", ".png"), width=1000, height=1000, type = "cairo-png")
  plotSmear(lrt, de.tags=deGenes)
  dev.off()
  srna_cnts_bin[[myF4]] <- lrt
}

# unique
for (myF4 in F4s_srnas) {
  message(myF4)
  # myF4 <-  F4s[7]
  # select only parents and F4s for time being:
  # should be same as above
  meta_unique_sub2 <- meta_unique_sub %>%
    dplyr::mutate(sub = (plant == myparent | plant == myF4))
  # DE analysis for unique merged sRNA bins
  group_unique <- subset(meta_unique_sub2, sub)$plant
  dgList_unique <- DGEList(counts=cnts_unique[,which(meta_unique_sub2$sub)],
                    group=group_unique,
                    remove.zeros = F,
                    genes=cnts_unique$ID
  )

  dgList_unique <- calcNormFactors(dgList_unique, method="TMM")
  design <- model.matrix(~group_unique)
  dgList_unique <- estimateDisp(dgList_unique, design)
  png(file=paste0(prefix, myF4, myparent,  "_unique", "_BCV", ".png"), width=1000, height=1000, type = "cairo-png")
  plotBCV(dgList_unique)
  dev.off()
  fit <- glmFit(dgList_unique,design)
  lrt <- glmLRT(fit, coef=2)
  deGenes <- decideTestsDGE(lrt, p=0.001)
  deGenes <- rownames(lrt)[as.logical(deGenes)]
  png(file=paste0(prefix, myF4, myparent, "_unique", "DE", ".png"), width=1000, height=1000, type = "cairo-png")
  plotSmear(lrt, de.tags=deGenes)
  dev.off()
  srna_cnts_unique_bin[[myF4]] <- lrt
}

map_int(srna_cnts_unique_bin, ~ sum(.x$table$PValue < 0.05))
#  P1512  P2561  P2562  P3611  P3612  P4041  P4042 
# 154417 164470 197502 181884 128041 176917 177524 



for (myF4 in F4s) {
  lrt <- srna_cnts_unique_bin[[myF4]]
  tags <- topTags(lrt, n = 1e5, p.value = 0.05)
  write_csv(as.data.frame(tags), path = paste0(prefix, mygenome, "_", myF4, "_table.csv"))
}
# }}}
save(srna_cnts_unique_bin, file = file.path(path_srnas, "analysis_bins" , paste0(prefix, mygenome,myparent, "_bins_raw_processed_unique.rdata")))

#------------------------------------ stats for dcl2 libraries {{{
# DESL for dcl2
meta_dcl2
#   File                                    index lib   plant   rep mapped generation lineage Run   include
#   <chr>                                   <dbl> <chr> <chr> <dbl> <chr>  <chr>      <chr>   <chr> <chr>  
# 1 wt-1_SolLyc_zeroMM_bowtie_unique.bam        3 wt-1  wt        1 SolLyc Parent     Parent  dcl2  yes    
# 2 wt-2_SolLyc_zeroMM_bowtie_unique.bam        4 wt-2  wt        2 SolLyc Parent     Parent  dcl2  yes    
# 3 dclab-1_SolLyc_zeroMM_bowtie_unique.bam     1 dcl-1 dcl2      1 SolLyc Parent     Parent  dcl2  yes    
# 4 dclab-2_SolLyc_zeroMM_bowtie_unique.bam     2 dcl-2 dcl2      2 SolLyc Parent     Parent  dcl2  yes    
dgList <- DGEList(counts=mycounts_dcl2[,1:4],
                  group=meta_dcl2$plant,
                  remove.zeros = F,
                  genes=mycounts_dcl2$ID
)
dgList <- calcNormFactors(dgList, method="TMM")
design <- model.matrix(~group)
dgList <- estimateDisp(dgList, design)
png(file=paste0(prefix, "dcl2", "_BCV", ".png"), width=1000, height=1000, type = "cairo-png")
plotBCV(dgList)
dev.off()
fit <- glmFit(dgList,design)
lrt_dcl2 <- glmLRT(fit, coef=2)
deGenes <- decideTestsDGE(lrt_dcl2, p=0.001)
deGenes <- rownames(lrt_dcl2)[as.logical(deGenes)]
png(file=paste0(prefix, "dcl2", "DE", ".png"), width=1000, height=1000, type = "cairo-png")
plotSmear(lrt_dcl2, de.tags=deGenes)
dev.off()

# load(file = paste0(prefix, mygenome, "_bins_raw_processed.rdata"))
# load(file = file.path(path_srnas, "analysis_bins" , paste0(prefix, mygenome, "_bins_raw_processed_unique.rdata")))

# }}}
#------------------------------------ stats for dcl2 libraries - WT vs WT infected  {{{
# DESL for dcl2
meta_dcl2
  # File           index lib   plant   rep mapped generation lineage Run   include
  # <chr>          <dbl> <chr> <chr> <dbl> <chr>  <chr>      <chr>   <chr> <chr>
# 1 wt-1_SolLyc_z…     3 wt-1  wt        1 SolLyc Parent     Parent  dcl2  yes
# 2 wt-2_SolLyc_z…     4 wt-2  wt        2 SolLyc Parent     Parent  dcl2  yes
# 3 wt-3_SolLyc_z…     4 wt-3  wt        3 SolLyc Parent     Parent  dcl2  yes
# 4 dclab-1_SolLy…     1 dcl-1 dcl2      1 SolLyc Parent     Parent  dcl2  yes
# 5 dclab-2_SolLy…     2 dcl-2 dcl2      2 SolLyc Parent     Parent  dcl2  yes
# 6 wt-1_tmv_SolL…     3 wt-1… wt_t…     1 SolLyc Parent     Parent  dcl2  yes
# 7 wt-2_tmv_SolL…     4 wt-2… wt_t…     2 SolLyc Parent     Parent  dcl2  yes
# 8 wt-3_tmv_SolL…     4 wt-3… wt_t…     3 SolLyc Parent     Parent  dcl2  yes
idx <- grep("wt", meta_dcl2$plant)
meta_sub <- meta_dcl2[idx,]
mycount_sub <- mycounts_dcl2[,idx]

dgList <- DGEList(counts=mycount_sub,
                  group=mycount_sub$plant,
                  remove.zeros = F,
                  genes=mycount_sub$ID
)
dgList <- calcNormFactors(dgList, method="TMM")
design <- model.matrix(~meta_sub$plant)
dgList <- estimateDisp(dgList, design)
png(file=paste0(prefix, "tmv", "_BCV", ".png"), width=1000, height=1000, type = "cairo-png")
plotBCV(dgList)
dev.off()
fit <- glmFit(dgList,design)
lrt_tmv <- glmLRT(fit, coef=2)
# deGenes <- decideTestsDGE(lrt_dcl2, p=0.001)
# deGenes <- rownames(lrt_dcl2)[as.logical(deGenes)]
png(file=paste0(prefix, "tmv", "DE", ".png"), width=1000, height=1000, type = "cairo-png")
plotSmear(lrt_dcl2, de.tags=deGenes)
dev.off()

# load(file = paste0(prefix, mygenome, "_bins_raw_processed.rdata"))
# load(file = file.path(path_srnas, "analysis_bins" , paste0(prefix, mygenome, "_bins_raw_processed_unique.rdata")))

# }}}
save(srna_cnts_bin, file = paste0(prefix, mygenome,myparent, "_bins_raw_processed.rdata"))
save(lrt_dcl2, mycounts_dcl2, file = paste0(prefix, mygenome, "_bins_raw_processed_dcl2.rdata"))
save(lrt_tmv, mycount_sub, file = paste0(prefix, mygenome, "_bins_raw_processed_tmv.rdata"))


# tmp
 load(file ="mirsl.rdata")
update.packages("tidyr")
library(tidyr)

 mirsl <- mirs %>%␠
   pivot_longer(cols = 2:23, names_to = "library", values_to = "counts")

   dplyr::select(miRNA, library, counts)

 mirsls <- mirsl %>%
   filter(!Parent_vs_P4042_PValue <= 0.05)

 mirslns <- mirsl %>%
   filter(!Parent_vs_P4042_PValue > 0.05)

 gg <- ggplot(mirsls, aes( y = counts, x = library)) +
   facet_grid(miRNA ~ ., scale = "free") +
   geom_bar(stat = 'identity', fill = 'grey') +
   gglayers2
 # geom_point()
 # theme(legend.position = "right",
 #       axis.text.x = element_text(angle = 300, hjust = 0),
 #       axis.ticks = element_blank()) +
 #       # scale_fill_manual(values = brewer.pal(4,"PRGn"))  +
 # scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6, digits = 2),
 #                    breaks = scales::pretty_breaks(n = 4))
 
 ggsave(paste0("miran_counts_counts_not_sigP4042vsParent", ".pdf"),
        gg,
        height = 70,
        width = 8,
        limitsize = F)

 gg <- ggplot(mirslns, aes( y = counts, x = library)) +
   facet_grid(miRNA ~ ., scale = "free") +
   geom_bar(stat = 'identity', fill = 'grey') +
   gglayers2

 ggsave(paste0("miran_counts_counts_non_sigP4042vsParent", ".pdf"),
        gg,
        height = 14,
        width = 8,
        limitsize = F)

 mirsls <- mirsl %>%
   filter(!Parent_vs_P4042_PValue < 0.05)
