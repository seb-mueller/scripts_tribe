
path_base <-  "/projects/TRIBE"
# laptop
path_base <-  "/home/sm934/workspace/tribe"
source(file.path(path_base, "scripts/project_settings.R"))
load(file=file.path(path_base, "current.workspace.Rdata"))
myF4 <- "P4042"
gtfilter <- "00"
genotype <- sym(paste0("genotype_", myF4))
mydf_filter4042 <- slydf %>%
  filter(!!genotype == gtfilter)


# mybarplot2(factor1, mycounts = count_list[[factor1]], mywidth = 5)


mydf_filter4042 <- mydf_filter4042 %>%
  select(-contains("P15")) %>%
  select(-contains("P36")) %>%
  select(-contains("PValue")) %>%
  select(-ends_with("class", ignore.case=F))

mydf_filter4042 <- mydf_filter4042 %>%
  mutate(hasExpression = srnas_P4042_logCPM < -2.79)

```

## Figures/facts DESL: P4042
```{r, DESLs}
load("/home/sm934/workspace/tribe/comparisons_using_big_table/sRNA_genomewide/final_figures_seb_code.rdata")

    # factor1 = "isPromotor"
    # factor2 = "chromatin_state"
    # factor3 = paste0(myprefix, myF4, infix, "_class")
    # context="CHH"
    myF4 ="P4042"
    infix = ""
    # gtfilter = "00"
    # myprefix = "srnas_"
    suffix = "_class2"
    # mydf = slydf
    myprefix <- "srnas_"
    suffix <- "_class2"
    factor2  <- paste0(myprefix, myF4, infix, suffix)


with(mydf_filter4042, table_prop(hasExpression)) #52% don't have any small RNAs (non loci): 3/4 of non_DESL(!)
mydf_filter4042  %>%
     dplyr::select(hasExpression, !!sym(factor2)) %>%
     table2()
#        srnas_P4042_class2
#          -DESL  +DESL inbetween non-DESL   <NA>
#   FALSE   3039  13843    401844   233613      0
#   TRUE       0      0         0   727538      0
#   <NA>       0      0         0        0      0


    mycounts <- mydf_filter4042 %>%
      count(hasExpression)

    gg <- ggplot(mycounts, aes(x = hasExpression, 
                           y = n)) +
          gglayers +
          geom_bar(stat = "identity") +
          geom_text(aes(label=n), vjust=-0.3) +
          scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3, digits = 2),
                             breaks = scales::pretty_breaks(n = 8))
    mytitle <- paste0("stats_P4042", "_sRNAloci", ".pdf")
    ggsave(mytitle, gg, height = 6, width = 2)

# ploting everything else
count_list <- list()

factors  <- c("annotation", "TEs_Order", "TEs_Order", "chromatin_state", 
              "TEs_Superfamily", "rnaseq_P4042_class2",
              "meth_P4042_CHH_class_stringent",
              "meth_P4042_CHG_class_stringent",
              "meth_P4042_CpG_class_stringent",
              "srnas_2124ratio_class",
              "TEs_rnaseq_P4042_class2"
)

for (factor1 in factors) {
count_list[[factor1]] <- mydf_filter4042 %>%
      filter(!!sym(factor1) != "none") %>%
      filter(!!sym(factor2) != "none") %>%
      count(!!sym(factor1), !!sym(factor2)) %>%
      group_by(!!sym(factor1)) %>%
      mutate(!!sym(paste0(factor1, "_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      group_by(!!sym(factor2)) %>%
      mutate(!!sym(paste0(factor2,"_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup()
}


for (factor1 in factors) {
mybarplot2(factor1, mycounts = count_list[[factor1]])
}


factor1 <- "annotation"
factor1 <- "TEs_Superfamily"
mybarplot2(factor1, mycounts = count_list[[factor1]], mywidth = 8)
factor1 <- "TEs_Order"
mybarplot2(factor1, mycounts = count_list[[factor1]], mywidth = 5)
factor1 <- "TEs_Class"
mybarplot2(factor1, mycounts = count_list[[factor1]], mywidth = 5)
# not much overlap DESL and TEs-DEGs

# more +DESLs
mydf_filter4042  %>% select(!!sym(factor2)) %>% table()
#     -DESL     +DESL inbetween  non-DESL 
#      3039     13843    401844    961151 
mydf_filter4042  %>% select(srnas_2124ratio_class) %>% table()
# high_24bp inbetween high_21bp 
#    137455   1210003     32419 
mydf_filter4042  %>%
     dplyr::select(srnas_2124ratio_class, !!sym(factor2)) %>%
     table2()
   #                      srnas_P4042_class2
   # srnas_2124ratio_class  -DESL  +DESL inbetween non-DESL   <NA>
   #             high_24bp   1066   3706     58299    74384      0
   #             inbetween   1548   5283    328585   874587      0
   #             high_21bp    425   4854     14960    12180      0
mydf_filter4042  %>%
     dplyr::select(srnas_2124ratio_class, !!sym(factor2)) %>%
     table_prop(mymargin=2)
# +DESL tend to be high in 21

mydf_filter4042  %>%
     dplyr::select(srnas_2124ratio_class, "TEs_Order")%>%
     table2()
mydf_filter4042  %>%
     dplyr::select(srnas_2124ratio_class, "TEs_Order")%>%
     table_prop(mymargin=2)
   #                      TEs_Order
   # srnas_2124ratio_class   DNA   LTR  none nonLTR Pararetrovirus   PHG Rolling   SSR    TE Unknown
   #             high_24bp 22.99  8.14  7.60  10.12           4.61 22.14    9.73 19.19 15.36   17.20
   #             inbetween 75.58 90.08 89.33  89.07          58.60 75.64   82.46 79.41 84.16   82.26
   #             high_21bp  1.43  1.78  3.08   0.80          36.79  2.21    7.80  1.40  0.48    0.54
# this could be a figure!!:
factor1="srnas_2124ratio_class"
factor2="TEs_Order"
mycounts <- slydf %>%
      count(!!sym(factor1), !!sym(factor2)) %>%
      group_by(!!sym(factor1)) %>%
      mutate(!!sym(paste0(factor1, "_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      group_by(!!sym(factor2)) %>%
      mutate(!!sym(paste0(factor2,"_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup()
mybarplot2(factor1,factor2, mycounts = mycounts, mywidth = 5)

factor2        = "srnas_P4042_class2"
mydf_filter4042  %>%
     dplyr::select(TEs_rnaseq_P4042_class2, !!sym(factor2))%>%
     table2()
   #                        srnas_P4042_class2
   # TEs_rnaseq_P4042_class2  -DESL  +DESL inbetween non-DESL   <NA>
   #               -DEG           1     11       122      214      0
   #               +DEG           4     62       679      916      0
   #               inbetween     36    584      4027     6270      0
   #               non-DEG      267   2524     69805   172177      0
   #               none        2731  10662    327211   781574      0
factor2="TEs_Order"
   # TEs_rnaseq_P4042_class2    DNA    LTR   none nonLTR Pararetrovirus    PHG Rolling    SSR     TE Unknown   <NA>
   #               -DEG           5    110      0      0             26      0       0      0      1      92    114
   #               +DEG         173    708      0     29             36      2       0      3     47      17    646
   #               inbetween    601   4060      0    138            274     27      87    107    116     456   5051
   #               non-DEG    10337  95733      0   2860           1472   1031     524   2438   9715   17187 103476
   #               none       50391 392129 484754  14279           8069   9611    1824  13039  45212  102829     41
mydf_filter4042  %>%
     dplyr::select(TEs_rnaseq_P4042_class2, !!sym(factor2))%>%
     table_prop(mymargin=2)
   # almost no overlap DESL and DETs

factor1="TEs_rnaseq_P4042_class2"
factor2="TEs_Order"
mycounts <- slydf %>%
      count(!!sym(factor1), !!sym(factor2)) %>%
      group_by(!!sym(factor1)) %>%
      mutate(!!sym(paste0(factor1, "_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      group_by(!!sym(factor2)) %>%
      mutate(!!sym(paste0(factor2,"_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup()
mybarplot2(factor1,factor2, mycounts = mycounts, mywidth = 5)

# save.image("/home/sm934/workspace/tribe/comparisons_using_big_table/sRNA_genomewide/final_figures_seb_code.rdata")
mybarplot(myF4, myprefix = "srnas_", factor1 = "chromatin_state", factor2 = "isPromotor", gtfilter = "00", mydf = slydf)
mybarplot(myF4, myprefix = "TEs_rnaseq_", infix = "", suffix = "_class2", factor1 = "srnas_P4042_class2", factor2 = "TEs_Order", gtfilter = "00", mydf = slydf, mywidth = 30)
mybarplot(myF4, myprefix = "TEs_rnaseq_", infix = "", suffix = "_class2", factor1 = "srnas_P4042_class2", factor2 = "TEs_Class", gtfilter = "00", mydf = slydf, mywidth = 30)
```
