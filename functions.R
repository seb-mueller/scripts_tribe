#------------------------------------ start here 
library(GenomicRanges)
library(datasets)
library(rtracklayer)
library(methylKit)
# library(genomation)
library(stringr)
library(Biostrings)
library(dplyr)
# library(topGO)
library(readr)
library(ggplot2)
library(purrr)
library(forcats)
library(viridis)
library(magrittr)
library(ggpubr)
# library(baySeq)
library(purrr)
try(library(rlang))

    gglayers_perc <- list(
                 scale_y_continuous(labels = scales::percent_format(acc =1),
                                    breaks = scales::pretty_breaks(n = 7)),
                 scale_fill_viridis_d(),
                 theme(legend.position = "right",
                       axis.ticks = element_blank(),
                       axis.text.x = element_text(angle = 300, hjust = 0)))

# mybarplot2(factor1, mycounts = count_list[[factor1]], mywidth = 5)
mybarplot2 <- function (factor1    = "TEs_Order",
                       factor2        = "srnas_P4042_class2",
                       mywidth        = 4,
                       mydf = mydf_filter4042,
                       ggplayers = gglayers,
                       mycounts = NULL) {

    gg <- ggplot(mycounts, aes(x = !!sym(factor2), 
                           y = !!sym(paste0(factor2, "_freq")), 
                           fill = fct_rev(!!sym(factor1)))) +
          gglayers +
          #           facet_grid(vars(), vars(!!sym(factor1))) +
          geom_bar(stat="identity", position="fill") + 
                 scale_fill_viridis_d() +
                 theme(legend.position = "right",
                       axis.ticks = element_blank(),
                       axis.text.x = element_text(angle = 300, hjust = 0))
    mytitle <- paste0(factor2, "_stats_rev_", factor1, ".pdf")
    ggsave(mytitle, gg, height = 6, width = mywidth)

    gg <- ggplot(mycounts, aes(x = !!sym(factor1), 
                           y = !!sym(paste0(factor2, "_freq")), 
                           fill = fct_rev(!!sym(factor2)))) +
#           gglayers +
          #           facet_grid(vars(), vars(!!sym(factor2))) +
          geom_bar(stat="identity", position="dodge") + 
                 scale_fill_viridis_d() +
                 theme(legend.position = "right",
                       axis.ticks = element_blank(),
                       axis.text.x = element_text(angle = 300, hjust = 0))
    mytitle <- paste0(factor2, "_stats_rev3_", factor1, ".pdf")
    ggsave(mytitle, gg, height = 6, width = mywidth+2)

    gg <- ggplot(mycounts, aes(x = !!sym(factor1), 
                           y = !!sym(paste0(factor1, "_freq")), 
                           fill = fct_rev(!!sym(factor2)))) +
          gglayers +
          #           facet_grid(vars(), vars(!!sym(factor1))) +
          geom_bar(stat="identity", position="fill") + 
                 scale_fill_viridis_d() +
                 theme(legend.position = "right",
                       axis.ticks = element_blank(),
                       axis.text.x = element_text(angle = 300, hjust = 0))
    mytitle <- paste0(factor2, "_stats_", factor1, ".pdf")
    ggsave(mytitle, gg, height = 6, width = mywidth)

    gg <- ggplot(mycounts, aes(x = !!sym(factor1), 
                           y = !!sym(paste0(factor1, "_freq")), 
                           fill = fct_rev(!!sym(factor2)))) +
#           gglayers +
          facet_grid(vars(), vars(!!sym(factor2))) +
          geom_bar(stat="identity", position="dodge") + 
                 scale_fill_viridis_d() +
                 theme(legend.position = "right",
                       axis.ticks = element_blank(),
                       axis.text.x = element_text(angle = 300, hjust = 0))
    mytitle <- paste0(factor2, "_stats_rev2_", factor1, ".pdf")
    ggsave(mytitle, gg, height = 6, width = mywidth+3)
}

mybarplot <- function (myF4           = "P4042",
                       gtfilter       = "00",
                       myprefix       = "srnas_",
                       infix          = "",
                       suffix         = "_class2",
                       factor1    = "TEs_Order",
                       factor2        = "isPromotor",
                       mywidth        = 12,
                       mydf = slydf) {
  require("ggpubr")
    # factor1 = "isPromotor"
    # factor2 = "chromatin_state"
    # factor3 = paste0(myprefix, myF4, infix, "_class")
    # context="CHH"
    # myF4 ="P4042"
    # infix = ""
    # gtfilter = "00"
    # myprefix = "srnas_"
    # suffix = "_class2"
    # mydf = slydf
    factor3  <- paste0(myprefix, myF4, infix, suffix)
    message(factor3)
    genotype <- sym(paste0("genotype_", myF4))

    mytitle <- paste0("stats_", factor1, "_", factor2, "_", 
                      factor3, "_filtered_", genotype, "-", gtfilter)
    gglayers <- list(
                 scale_y_continuous(labels = scales::percent,
                                    breaks = scales::pretty_breaks(n = 8)),
                 scale_fill_viridis_d(),
                 theme(legend.position = "right",
                       axis.ticks = element_blank(),
                       axis.text.x = element_text(angle = 300, hjust = 0)))

    mydf_filter <- mydf %>%
      filter(!!genotype == gtfilter)

    mycounts <- mydf_filter %>%
      count(!!sym(factor1))
    gga <- ggplot(mycounts, aes(x = !!sym(factor1), 
                           y = n)) +
          gglayers +
          geom_bar(stat = "identity") +
          geom_text(aes(label=n), vjust=-0.3) +
          scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3, digits = 2),
                             breaks = scales::pretty_breaks(n = 8))
    mycounts <- mydf_filter %>%
      count(!!sym(factor2))
    ggb <- ggplot(mycounts, aes(x = !!sym(factor2), 
                           y = n)) +
          gglayers +
          geom_bar(stat = "identity") +
          geom_text(aes(label=n), vjust=-0.3) +
          scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3, digits = 2),
                             breaks = scales::pretty_breaks(n = 8))
    mycounts <- mydf_filter %>%
      count(!!sym(factor3))
    ggc <- ggplot(mycounts, aes(x = !!sym(factor3), 
                           y = n)) +
          gglayers +
          geom_bar(stat = "identity") +
          geom_text(aes(label=n), vjust=-0.3) +
          scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3, digits = 2),
                             breaks = scales::pretty_breaks(n = 8))

    mycounts <- mydf_filter %>%
      # group_modify(~as.data.frame((table(!!sym(myplantsclass)))))
      filter(!!sym(factor1) != "none") %>%
      filter(!!sym(factor2) != "none") %>%
      filter(!!sym(factor3) != "none") %>%
      count(!!sym(factor1), !!(sym(factor2)), !!(sym(factor3))) %>%
      group_by(!!sym(factor1)) %>%
      mutate(!!sym(paste0(factor1, "_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      group_by(!!sym(factor2)) %>%
      mutate(!!sym(paste0(factor2,"_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup() %>%
      group_by(!!sym(factor3)) %>%
      mutate(!!sym(paste0(factor3, "_freq")) := round(100 * n / sum(n),2)) %>%
      ungroup()

    gg0 <- ggplot(mycounts, aes(x = !!sym(factor1), 
                           y = n)) +
          facet_grid(vars(), vars(!!sym(factor2), fct_rev(!!sym(factor3)))) +
          gglayers +
          geom_bar(stat = "identity") +
          scale_y_continuous(labels = scales::unit_format(unit = "k", scale = 1e-3, digits = 2),
                             breaks = scales::pretty_breaks(n = 8)) +
                     labs(title = mytitle)

    # mytitle <- paste0("barplot_x-", factor1, "_y-", factor2, "_fill-", 
                      # factor3, "_filtered_", genotype, "-", gtfilter, ".pdf")
    # message(mytitle)
    gg1 <- ggplot(mycounts, aes(x = !!sym(factor1), 
                           y = !!sym(paste0(factor1, "_freq")), 
                           fill = fct_rev(!!sym(factor3)))) +
          facet_grid(vars(), vars(!!sym(factor3), !!sym(factor2))) +
          geom_bar(stat="identity", position="dodge") + 
                 scale_fill_viridis_d() +
                 theme(legend.position = "right",
                       axis.ticks = element_blank(),
                       axis.text.x = element_text(angle = 300, hjust = 0))
          # gglayers
    # ggsave(mytitle, gg, height = 6, width = 6)

    # mytitle <- paste0("barplot_x-", factor3, "_y-", factor2, "_fill-", 
    #                   factor1, "_filtered_", genotype, "-", gtfilter, ".pdf")
    # message(mytitle)
     # browser()

    gg2 <- ggplot(mycounts, aes(x = !!sym(factor3), 
                           fill = fct_rev(!!sym(paste0(factor1))), 
                           y = !!sym(paste0(factor3, "_freq")))) +
          facet_grid(vars(), vars(!!sym(factor1), !!sym(factor2))) +
          geom_bar(stat="identity", position="dodge")  +
                 scale_fill_viridis_d() +
                 theme(legend.position = "right",
                       axis.ticks = element_blank(),
                       axis.text.x = element_text(angle = 300, hjust = 0))
          # gglayers
    # ggsave(mytitle, gg, height = 6, width = 6)

    # mytitle <- paste0("barplot_x-", factor3, "_y-", factor1, "_fill-", 
    #                   factor2, "\n_filtered_", genotype, "-", gtfilter, ".pdf")
    # message(mytitle)
    gg3 <- ggplot(mycounts, aes(x = !!sym(factor3), 
                           y = !!sym(paste0(factor3, "_freq")), 
                           fill = fct_rev(!!sym(factor2)))) +
          facet_grid(vars(), vars(!!sym(factor1))) +
          geom_bar(stat="identity", position="fill") + 
          gglayers

    # ggsave(mytitle, gg, height = 6, width = 6)
        ggall <- ggarrange(gg0, 
                   ggarrange(ggarrange(gga, ggb, ggc, ncol = 3), gg1), 
                   ggarrange(gg2, gg3), nrow = 3)
                    # labels = c("A", "B", "C", "D"),
                    # ncol = 2, nrow = 2)
    mytitle <- paste0(factor3, "_", factor1, "_", 
                      factor2, "_", genotype, "-", gtfilter)
    ggsave(paste0(mytitle, ".pdf"), ggall, height = 15, width = mywidth)
    write_tsv(mycounts, path = paste0(mytitle, ".tsv"))
    print(paste0(mytitle, ".pdf"))
    return(ggall)
}
