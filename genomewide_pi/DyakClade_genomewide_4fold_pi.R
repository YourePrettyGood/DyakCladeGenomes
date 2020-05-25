#Load libraries to use:
library(tidyverse)

#Set working directory:
setwd('/home/pfreilly/Documents/Manuscripts/DyakCladeGenomes/DyakCladeGenomes/genomewide_pi/')

#This script expects that `make -f calculate_pi.mk all` has already been run

#General variables:
chrom_arms <- c("X", "2L", "2R", "3L", "3R")
chromosomes <- c("X", "2", "3")
species <- c("Dsan", "Dyak", "Dtei")
gap_length <- 500L

#Plotting details:
linesize <- 0.5
fontsize <- 12
keysize <- 30
ytickfontsize <- 8
tmargintrim <- 10
bmargintrim <- 2
width_cm <- 16.0
width_in <- width_cm / 2.54
height_cm <- 6.0
height_in <- height_cm / 2.54

#Load in the 4-fold pi estimates for 100 kb windows for each species:
pi_col_classes <- c("character", "integer", "numeric", "numeric")
pi_col_names <- c("Scaffold", "WindowStart", "Pi4f", "UsableFraction")
Dsan_pi <- read.table(gzfile('Dsan_4fold_pi_w100000.tsv.gz'),
                      colClasses=pi_col_classes, col.names=pi_col_names)
Dtei_pi <- read.table(gzfile('Dtei_4fold_pi_w100000.tsv.gz'),
                      colClasses=pi_col_classes, col.names=pi_col_names)
Dyak_pi <- read.table(gzfile('Dyak_4fold_pi_w100000.tsv.gz'),
                      colClasses=pi_col_classes, col.names=pi_col_names)
Dsan_pi$Species <- "Dsan"
Dtei_pi$Species <- "Dtei"
Dyak_pi$Species <- "Dyak"
joint_pi_df <- rbind(Dsan_pi, Dtei_pi, Dyak_pi)

#Load in the FASTA indexes for the reference genomes so that we can
# consolidate the L and R arms into full chromosomes:
fai_col_classes <- c("character", rep("integer", 4))
fai_col_names <- c("Scaffold", "Length", "Offset", "Wrap", "Wrap1")
Dsan_fai <- read.table('FAIs/Dsan.fasta.fai',
                       colClasses=fai_col_classes, col.names=fai_col_names)
Dtei_fai <- read.table('FAIs/Dtei.fasta.fai',
                       colClasses=fai_col_classes, col.names=fai_col_names)
Dyak_fai <- read.table('FAIs/Dyak.fasta.fai',
                       colClasses=fai_col_classes, col.names=fai_col_names)

#Now we have to manipulate the pi data.frame a bit to consolidate L and R
# arms into full chromosomes, with a 500 bp gap between arms, leaving all
# other sequences alone, and converting coordinates from bp to Mbp:
joint_pi_df <- joint_pi_df %>%
   mutate(Chromosome=case_when(Scaffold == "X" ~ "X",
                               Scaffold %in% c("2L", "2R") ~ "2",
                               Scaffold %in% c("3L", "3R") ~ "3",
                               ! Scaffold %in% chrom_arms ~ Scaffold),
          WindowPos=case_when(Scaffold == "2R" & Species == "Dsan" ~ WindowStart + subset(Dsan_fai, Scaffold == "2L")$Length + gap_length,
                              Scaffold == "2R" & Species == "Dtei" ~ WindowStart + subset(Dtei_fai, Scaffold == "2L")$Length + gap_length,
                              Scaffold == "2R" & Species == "Dyak" ~ WindowStart + subset(Dyak_fai, Scaffold == "2L")$Length + gap_length,
                              Scaffold == "3R" & Species == "Dsan" ~ WindowStart + subset(Dsan_fai, Scaffold == "3L")$Length + gap_length,
                              Scaffold == "3R" & Species == "Dtei" ~ WindowStart + subset(Dtei_fai, Scaffold == "3L")$Length + gap_length,
                              Scaffold == "3R" & Species == "Dyak" ~ WindowStart + subset(Dyak_fai, Scaffold == "3L")$Length + gap_length,
                              ! Scaffold %in% c("2R", "3R") ~ WindowStart)/1e6)

#To a certain extent, this hides the variation in pericentromeric assembly,
# e.g. Dsan has a longer 2L assembly than Dyak and especially Dtei, but
# Dtei has a longer 2R assembly than the other two.
#However, this gets more to the point of the figure: Patterns of diversity
# across the chromosomes between the three species.

#Plot the figure:
joint_pi_df %>%
   filter(Chromosome %in% chromosomes) %>%
   mutate(Chromosome=factor(Chromosome, levels=chromosomes, labels=chromosomes),
          Species=factor(Species, levels=species, labels=species),
          Pi4f=Pi4f*100) %>%
   ggplot(aes(x=WindowPos, y=Pi4f, colour=Species)) +
      geom_line(size=linesize) +
      facet_grid(Species ~ Chromosome, scales='free_x', space='free_x') +
      scale_x_continuous(breaks=seq(0, 60, 10)) +
      scale_y_continuous(breaks=seq(0.0, 6.0, 3.0)) +
      labs(x="Position (Mbp)", y=expression(pi["4f"]~"(%)")) +
      theme_bw() +
      theme(text=element_text(size=fontsize),
            legend.title=element_text(size=fontsize, face="bold"),
            axis.title=element_text(size=fontsize, face="bold"),
            plot.title=element_text(size=fontsize, face="bold"),
            axis.text.y=element_text(size=ytickfontsize),
            legend.position="top",
            legend.justification="left",
            legend.key.width=unit(keysize, 'pt'),
            legend.box.margin=margin(t=-tmargintrim, r=0, b=-bmargintrim, l=0, unit='pt'),
            legend.box.spacing=unit(2, 'pt'),
            strip.text.y=element_blank(),
            panel.spacing=unit(0, 'cm'))

ggsave('DyakClade_genomewide_4fold_pi_w100kb.pdf', width=width_cm, height=height_cm, units="cm", dpi=500)

combined_depth_df <- read.table(gzfile('/home/pfreilly/Documents/Manuscripts/DyakCladeGenomes/DyakCladeGenomes/genomewide_pi/combined_4fold_depth_w100kb.tsv.gz'), colClasses=c("character", "integer", "character", "character", "character", "numeric", "numeric"), col.names=c("Scaffold", "WindowStart", "Statistic", "Species", "Line", "Depth", "UsableFraction"))
#Now we have to manipulate the depth data.frame a bit to consolidate L and R
# arms into full chromosomes, with a 500 bp gap between arms, leaving all
# other sequences alone, and converting coordinates from bp to Mbp:
combined_depth_df <- combined_depth_df %>%
   mutate(Chromosome=case_when(Scaffold == "X" ~ "X",
                               Scaffold %in% c("2L", "2R") ~ "2",
                               Scaffold %in% c("3L", "3R") ~ "3",
                               ! Scaffold %in% chrom_arms ~ Scaffold),
          WindowPos=case_when(Scaffold == "2R" & Species == "Dsan" ~ WindowStart + subset(Dsan_fai, Scaffold == "2L")$Length + gap_length,
                              Scaffold == "2R" & Species == "Dtei" ~ WindowStart + subset(Dtei_fai, Scaffold == "2L")$Length + gap_length,
                              Scaffold == "2R" & Species == "Dyak" ~ WindowStart + subset(Dyak_fai, Scaffold == "2L")$Length + gap_length,
                              Scaffold == "3R" & Species == "Dsan" ~ WindowStart + subset(Dsan_fai, Scaffold == "3L")$Length + gap_length,
                              Scaffold == "3R" & Species == "Dtei" ~ WindowStart + subset(Dtei_fai, Scaffold == "3L")$Length + gap_length,
                              Scaffold == "3R" & Species == "Dyak" ~ WindowStart + subset(Dyak_fai, Scaffold == "3L")$Length + gap_length,
                              ! Scaffold %in% c("2R", "3R") ~ WindowStart)/1e6)

joint_pi_df %>%
   full_join(combined_depth_df,
             by=c("Chromosome", "WindowPos", "Species", "Scaffold", "WindowStart")) %>%
   filter(Chromosome %in% chromosomes) %>%
   mutate(Chromosome=factor(Chromosome, levels=chromosomes, labels=chromosomes),
          Species=factor(Species, levels=species, labels=species),
          Pi4f=Pi4f*100) %>%
   ggplot(aes(x=WindowPos, y=Pi4f, colour=Species)) +
      geom_line(size=linesize) +
      geom_quantile(aes(y=Depth/10), linetype=4, method="rqss", lambda=1) +
      facet_grid(Species ~ Chromosome, scales='free_x', space='free_x') +
      scale_x_continuous(breaks=seq(0, 60, 10)) +
      scale_y_continuous(breaks=seq(0.0, 9.0, 3.0), limits=c(0.0, 10.0)) +
      labs(x="Position (Mbp)", y=expression(pi["4f"]~"(%)"~"and"~"depth/10")) +
      theme_bw() +
      theme(text=element_text(size=fontsize),
            legend.title=element_text(size=fontsize, face="bold"),
            axis.title=element_text(size=fontsize, face="bold"),
            plot.title=element_text(size=fontsize, face="bold"),
            axis.text.y=element_text(size=ytickfontsize),
            legend.position="top",
            legend.justification="left",
            legend.key.width=unit(keysize, 'pt'),
            legend.box.margin=margin(t=-tmargintrim, r=0, b=-bmargintrim, l=0, unit='pt'),
            legend.box.spacing=unit(2, 'pt'),
            strip.text.y=element_blank(),
            panel.spacing=unit(0, 'cm'))

ggsave('DyakClade_genomewide_4fold_pi_wDepthQuantiles_w100kb.pdf', width=width_cm, height=height_cm*2, units="cm", dpi=500)

#Quick reset:
rm(chrom_arms, chromosomes, species, gap_length)
rm(linesize, fontsize, keysize, ytickfontsize, tmargintrim, bmargintrim, width_in, width_cm, height_in, height_cm)
rm(pi_col_classes, pi_col_names, Dsan_pi, Dtei_pi, Dyak_pi, joint_pi_df)
rm(fai_col_classes, fai_col_names, Dsan_fai, Dtei_fai, Dyak_fai)
rm(combined_depth_df)
