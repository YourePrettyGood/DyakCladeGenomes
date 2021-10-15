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
#BUSCO bar height:
busco_bar_height <- 0.75
#Plot title:
busco_title <- "BUSCO Assembly Completeness Comparison"
#Colour palette:
#colorblind_palette <- c("#009E73", "#E79F00", "#9AD0F3", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#000000")
busco_colors <- c("#56B4E9", "#3492C7", "#F0E442", "#F04442")
#Plot dimensions (in cm):
busco_dims <- c(16, 12)
#Normally, BUSCO would specify the font family, font size ratio, and labelsize:
#busco_fontfamily <- "sans"
#busco_fontsizeratio <- 1
#busco_labsize <- 1
#if (length(unique(busco_df$Species)) > 10) {
#   busco_labsize <- 0.66
#}

#Labels and level order for BUSCO classes:
busco_levels <- c("S", "D", "F", "M")
busco_labels <- c("Complete (C) and single-copy (S)",
                  "Complete (C) and duplicated (D)",
                  "Fragmented (F)",
                  "Missing (M)")

#BUSCO does a silly thing and hard-codes the columns of the data.frame using
# generate_plot.py, so we avoid that and save the lines of a table using
# calls to BUSCOsummary.awk
#Read in the BUSCO results compiled using BUSCOsummary.awk:
busco_df <- read.table(input_file, header=TRUE)

#Set the assembly order for plotting (has to be reversed to plot the way
# you'd expect):
busco_asms <- rev(asm_names)
asm_labels <- gsub("_", " ", busco_asms, fixed=TRUE)

#Generate the BUSCO text strings to be inlaid in the bars:
busco_strings <- busco_df %>%
   mutate(Species=factor(Species,
                         levels=busco_asms,
                         labels=asm_labels)) %>%
   dplyr::select(-Percentage) %>%
   pivot_wider(id_cols=Species,
               names_from=Category,
               values_from=Count) %>%
   group_by(Species) %>%
   summarize(BUSCO=paste0("C:", C,
                          " [S:", S,
                          ", D:", D,
                          "], F:", F,
                          ", M:", M,
                          ", n:", N))

#Generate the BUSCO plot:
busco_plot <- busco_df %>%
   mutate(Species=factor(Species,
                         levels=busco_asms,
                         labels=asm_labels)) %>%
   filter(! Category %in% c("C", "N")) %>%
   mutate(Category=factor(Category,
                          levels=busco_levels,
                          labels=busco_labels)) %>%
   ggplot(aes(x=Species, y=Percentage, fill=Category)) +
      geom_col(width=busco_bar_height) +
      coord_flip() +
      annotate("text",
               label=busco_strings$BUSCO,
               x=seq(1, nrow(busco_strings)),
               y=99,
               size=4,
               colour="black",
               hjust=1) +
      theme_bw() +
      scale_y_continuous(breaks=seq(0, 100, by=20)) +
      scale_fill_manual(values=busco_colors, labels=busco_labels) +
      labs(x="",
           y="% BUSCOs",
           title=busco_title) +
      theme(legend.position="top",
            legend.title=element_blank(),
            panel.background=element_rect(color="#FFFFFF", fill="white"),
            panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            axis.ticks.y=element_blank(),
            axis.line = element_line(colour = "black"),
            text=element_text(size=12),
            legend.text=element_text(size=8),
            axis.title.y=element_text(vjust=1,
                                      face="bold",
                                      size=12,
                                      margin=margin(10,20,10,0)),
            axis.title.x=element_text(vjust=1,
                                      face="bold",
                                      size=12,
                                      margin=margin(10,20,10,0)),
            axis.text.x=element_text(size=11),
            axis.text.y=element_text(size=11),
            plot.title=element_text(face="bold")) +
      guides(fill=guide_legend(override.aes=list(colour=NULL),
                               nrow=2,
                               byrow=TRUE))

ggsave(plot_file,
       busco_plot,
       width=busco_dims[1],
       height=busco_dims[2],
       units="cm",
       dpi=500)
