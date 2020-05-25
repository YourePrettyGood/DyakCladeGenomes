#!/usr/bin/env Rscript

options <- commandArgs(trailingOnly=TRUE)

check_package <- function(pkg_name) {
   if (!require(pkg_name, character.only=TRUE)) {
      install.packages(pkg_name, dependencies=TRUE, type="source", repos="http://cran.us.r-project.org")
   }
   library(pkg_name, character.only=TRUE)
}

check_package("tidyverse")

#Read in arguments:
plot_file <- options[1]
input_file <- options[2]
asm_names <- options[-c(1,2)]

#Plotting details:
colorblind_palette <- c("#009E73", "#E79F00", "#9AD0F3", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#000000")

all_lens <- read.table(input_file, colClasses=c("character", "integer", "integer"), col.names=c("Assembly", "Index", "Length"), header=FALSE)

Assembly_levels <- asm_names
Assembly_labels <- gsub("_", " ", asm_names, fixed=TRUE)
#Assembly_levels <- c("Dmel_ISO1_FBr6.26_scafs", "Dmel_ISO1_FBr6.26", "Dsan_STOCAGO1482", "Dtei_GT53w", "Dyak_NY73PB", "Dyak_Tai18E2", "Dyak_Tai18E2_FBr1.05")
#Assembly_labels <- c("Ideal Dmel ISO1", "Dmel ISO1", "Dsan STO CAGO 1482", "Dtei GT53w", "Dyak NY73PB", "Dyak Tai18E2", "Dyak Tai18E2 FlyBase")
contiguity_plot <- all_lens %>% mutate(Assembly=factor(Assembly, levels=Assembly_levels, labels=Assembly_labels)) %>%
  group_by(Assembly) %>%
  mutate(CumulativeLength=cumsum(Length)) %>%
  ggplot(aes(x=Index, y=CumulativeLength, colour=Assembly)) +
    geom_line(size=1) +
    scale_x_log10() +
    scale_colour_manual(values=colorblind_palette) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
      panel.border=element_rect(colour="black", size=1.5),
      aspect.ratio=0.8,
      text=element_text(size=12),
      legend.text=element_text(size=8),
      legend.title=element_text(face="bold", size=12),
      axis.title.y=element_text(vjust=1, face="bold", size=12, margin=margin(10,20,10,0)),
      axis.title.x=element_text(vjust=1, margin=margin(10,20,10,0), face="bold", size=12),
      axis.text.x=element_text(size=11),
      axis.text.y=element_text(size=11),
      plot.title=element_text(face="bold")) +
  labs(x=expression(paste("Contig Rank (", log[10], " scale)")),
    y="Cumulative Assembly Length (bp)",
    title="Assembly Contiguity")
ggsave(plot_file, plot=contiguity_plot, dpi=500, units="cm", width=16.0, height=12.0)
