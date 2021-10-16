#!/usr/bin/env Rscript

#Get the arguments from the command line:
options <- commandArgs(trailingOnly=TRUE)
offsets_file <- options[1]
step_size <- options[2]
window_sizes <- options[-seq(1,2)]

#Load libraries to use:
library(tidyverse)

#Set working directory:
setwd('/home/pfreilly/Documents/Manuscripts/DyakCladeGenomes/DyakCladeGenomes/HKA/')

#Plotting details:
linesize <- 0.5
fontsize <- 12
keysize <- 30
ytickfontsize <- 8
tmargintrim <- 10
bmargintrim <- 2
width_cm <- 16.0
width_in <- width_cm / 2.54
height_cm <- 12.0
height_in <- height_cm / 2.54

#Analysis details:
#For now we hard-code these:
chrom_arms <- c("X", "2L", "2R", "3L", "3R")
species <- c("Dsan", "Dyak", "Dtei", "Dmel", "Dsim")
#window_sizes <- c(1000000)
#step_size <- 200000
avg_types <- c("weighted", "naive")
site_class <- "4fold"

#Functions for automating read-in:
#We use a file of offsets to "align" the arms of the different species,
# since there's variation in contiguity between assemblies:
read.offsets <- function(fn) {
   read.table(fn, colClasses=c("character", "character", "character", "numeric"), header=TRUE) %>%
      mutate(Scaffold=str_replace(Scaffold, "^Scf_", "")) #Just for Dsim, remove the Scf_ prefix for scaffolds
}

#We use this BED to transform between window IDs and window start positions:
read.window.bed <- function(species, window_size, step_size) {
   fn <- paste0("BEDs/", species, "_w", window_size, "_s", step_size, "_windows.bed");
   read.table(fn, colClasses=c("character", "integer", "integer"),
              col.names=c("Scaffold", "BEDWindowStart", "WindowEnd"),
              header=FALSE) %>%
      mutate(WindowStart=BEDWindowStart+1,
             Species=species) %>%
      select(-BEDWindowStart) %>%
      mutate(Scaffold=str_replace(Scaffold, "^Scf_", "")) #Just for Dsim, remove the Scf_ prefix for scaffolds
}

#Read in the table with observed and bootstrap values per window,
# and be sure to annotate which file it came from for collation:
read.bootstrap.table <- function(species, window_size, step_size, site_class, statistic, avg_type) {
   fn <- paste0("stats/", species, "_w", window_size, "_s", step_size, "_", statistic, "_", site_class, "_", avg_type, "Avg_obs_wBootstraps.tsv.gz");
   read.table(gzfile(fn), header=1, stringsAsFactors=FALSE) %>%
      mutate(Species=!!species,
             Statistic=!!statistic,
             AvgType=!!avg_type,
             WindowSize=!!window_size)
}

#Rename the value column and remove excess columns ahead of the join:
trim.excess.columns <- function(df, statistic) {
   omit_cols <- c("Numerator", "Denominator", "NumOGs", statistic)
   df %>% mutate(Value=!!sym(statistic)) %>%
      select(-!!omit_cols)
}

#This function performs the transformations of positions for plotting.
# These transformations include:
#  1) Converting from Window IDs to window start positions (1-based)
#  2) Concatenating the L and R arms of each chromosome
#  3) Scaling the positions so the numbers are easy to plot (Default: bp to Mbp)
#The input data.frames must contain a single species -- collation happens
# later.
#The function is custom to Drosophila with the same overall karyotype as
# D. melanogaster.
transform.coordinates <- function(df, bed, offsets, chrom_arms, bp_scale_factor=1e6) {
   stopifnot(bp_scale_factor > 0);
   df %>% mutate(WindowStart=bed[WindowID, "WindowStart"],
                 Scaffold=bed[WindowID, "Scaffold"]) %>%
      left_join(offsets, by=c("Scaffold", "Species")) %>%
      mutate(WindowPos=case_when(Scaffold %in% chrom_arms ~ WindowStart + Offset,
                                 ! Scaffold %in% chrom_arms ~ WindowStart)) %>%
      mutate(WindowPos=WindowPos/bp_scale_factor) %>%
      select(-WindowStart)
}

#This function merges the Pi and Dxy tibbles and calculates Pi/Dxy
# grouped on AvgType, Species, WindowID, Chromosome, and WindowPos.
#Expects input to be tall with respect to statistics and their values,
# i.e. statistic name in column "Statistic" and value in column "Value".
pi.over.Dxy <- function(pi, dxy, site_class) {
   bind_rows(pi, dxy) %>% spread(Statistic, Value) %>%
      mutate(PiOverDxy=Pi/Dxy) %>%
      gather("Statistic", "Value", Pi, Dxy, PiOverDxy)
}

#This function identifies the two-tailed bootstrap confidence interval
# for each window using the quantile method.
#Note that the input needs to be tall with respect to statistics and
# their values, in columns named "Statistic" and "Value, respectively.
bootstrap.CI.quantile <- function(df, coverage) {
   stopifnot(coverage > 0, coverage < 1);
   upper_quantile <- (1+coverage)/2;
   lower_quantile <- (1-coverage)/2;
   bootstrap_df <- df %>% filter(ID > 0) %>%
      group_by(Statistic, AvgType, WindowSize, Species, WindowID, Scaffold, WindowPos) %>%
      summarize(lowerCI=quantile(Value, probs=lower_quantile, na.rm=TRUE),
                upperCI=quantile(Value, probs=upper_quantile, na.rm=TRUE))
   df %>% filter(ID == 0) %>%
      select(-ID) %>%
      full_join(bootstrap_df,
                by=c("Statistic", "AvgType", "WindowSize", "Species",
                     "WindowID", "Scaffold", "WindowPos"))
}

#This function fills in missing windows with NAs so that we don't get
# funky lines drawn in NA regions:
add.missing.windows <- function(df, bed, offsets, chrom_arms, bp_scale_factor=1e6) {
   all_combinations <- expand_grid(WindowID=seq(1, nrow(bed)),
                                   Species=unique(df$Species),
                                   AvgType=unique(df$AvgType),
                                   WindowSize=unique(df$WindowSize),
                                   Statistic=unique(df$Statistic)) %>%
      transform.coordinates(bed=bed, offsets=offsets, chrom_arms=chrom_arms, bp_scale_factor=bp_scale_factor);
   df %>% right_join(all_combinations)
}

stats_wCIs_list <- list()
offsets <- read.offsets(offsets_file)
for (avg_type in avg_types) {
   for (window_size in window_sizes) {
      for (spp in species) {
         cat(paste("Processing", window_size, avg_type, "averages for species", spp, "\n"))
         bed <- read.window.bed(species=spp,
                                window_size=window_size,
                                step_size=step_size)
         pi <- read.bootstrap.table(species=spp,
                                    window_size=window_size,
                                    step_size=step_size,
                                    site_class=site_class,
                                    statistic="Pi",
                                    avg_type=avg_type) %>%
            trim.excess.columns("Averagepi_4Fold") %>%
            transform.coordinates(bed=bed,
                                  offsets=offsets,
                                  chrom_arms=chrom_arms)
         dxy <- read.bootstrap.table(species=spp,
                                     window_size=window_size,
                                     step_size=step_size,
                                     site_class=site_class,
                                     statistic="Dxy",
                                     avg_type=avg_type) %>%
            trim.excess.columns("AverageDxy_4Fold") %>%
            transform.coordinates(bed=bed,
                                  offsets=offsets,
                                  chrom_arms=chrom_arms)
         all_stats <- pi.over.Dxy(pi=pi, dxy=dxy, site_class=site_class)
         stats_wCIs_list[[length(stats_wCIs_list)+1]] <- all_stats %>% 
            bootstrap.CI.quantile(coverage=0.95) %>%
            add.missing.windows(bed=bed, offsets=offsets, chrom_arms=chrom_arms)
      }
   }
}

#Concatenate the results from the different average types, window sizes, and
# species:
stats_wCIs <- bind_rows(stats_wCIs_list)

#Save the point estimates with bootstrap CI bounds as an Rdata file
# for Kevin to use:
save(stats_wCIs, file='HKA_4fold_stats_wCIs.Rdata',
     compress=TRUE, compression_level=9)

#Now we plot pi and Dxy on the two y axes with CIs:
for (avg_type in avg_types) {
   for (window_size in window_sizes) {
      #Set the ymax for pi:
      ymax <- 0.06
      #and the Dxy scale factor (effectively, Dxy's ymax will be this *times* ymax):
      dxy_scale_factor <- 5
      pi_and_dxy_plot <- stats_wCIs %>%
         filter(Scaffold %in% chrom_arms,
                Statistic %in% c("Pi", "Dxy"),
                AvgType == avg_type,
                WindowSize == window_size) %>%
         mutate(Value=case_when(Statistic == "Dxy" ~ Value/!!dxy_scale_factor,
                                Statistic != "Dxy" ~ Value),
                lowerCI=case_when(Statistic == "Dxy" ~ lowerCI/!!dxy_scale_factor,
                                  Statistic != "Dxy" ~ lowerCI),
                upperCI=case_when(Statistic == "Dxy" ~ upperCI/!!dxy_scale_factor,
                                  Statistic != "Dxy" ~ upperCI),
                Scaffold=factor(Scaffold, levels=chrom_arms),
                Species=factor(Species, levels=species)) %>%
         ggplot(aes(x=WindowPos, y=Value,
                    colour=Statistic, fill=Statistic,
                    ymin=lowerCI, ymax=upperCI)) +
            geom_line(size=0.25) +
            geom_ribbon(linetype=0, alpha=0.5) +
            theme_bw() +
            labs(x="Position (Mbp)",
                 title=paste0("window=", as.numeric(window_size)/1e3, " kbp, step=", as.numeric(step_size)/1e3, " kbp")) +
            scale_x_continuous(breaks=seq(0, 60, by=10)) +
            scale_y_continuous(limits=c(0, ymax),
                               name=expression(pi),
                               sec.axis=sec_axis(~ . * dxy_scale_factor,
                                                 name=expression(D[xy]))) +
            scale_colour_manual(values=c("blue", "black")) +
            scale_fill_manual(values=c("blue", "black")) +
            guides(colour=FALSE, fill=FALSE) +
            facet_grid(Species ~ Scaffold,
                       scales="free_x",
                       space="free_x") +
            theme(axis.title=element_text(face="bold"),
                  panel.spacing.x=unit(0, "cm"),
                  text=element_text(size=fontsize),
                  axis.text=element_text(size=fontsize-2),
                  axis.title.y.left=element_text(colour="black"),
                  axis.title.y.right=element_text(colour="blue"),
                  axis.ticks.y.right=element_line(colour="blue"),
                  axis.text.y.right=element_text(colour="blue"))
#      print(pi_and_dxy_plot)
      ggsave(paste0('HKA_4fold_', avg_type, 'Avg_w', window_size, '_pi_and_Dxy.pdf'),
             plot=pi_and_dxy_plot,
             width=width_cm,
             height=height_cm,
             units="cm",
             dpi=500)

      #Make sure to set the ymax appropriately for pi/Dxy:
      ymax <- 0.35
      #Plot the HKA-like ratio with CIs:
      pi_over_dxy_plot <- stats_wCIs %>%
         filter(Scaffold %in% chrom_arms,
                Statistic %in% c("PiOverDxy"),
                AvgType == avg_type,
                WindowSize == window_size) %>%
         mutate(Scaffold=factor(Scaffold, levels=chrom_arms),
                Species=factor(Species, levels=species)) %>%
         ggplot(aes(x=WindowPos, y=Value,
                    colour=Statistic, fill=Statistic,
                    ymin=lowerCI, ymax=upperCI)) +
            geom_line(size=0.25) +
            geom_ribbon(linetype=0, alpha=0.5) +
            theme_bw() +
            labs(x="Position (Mbp)",
                 title=paste0("window=", as.numeric(window_size)/1e3, " kbp, step=", as.numeric(step_size)/1e3, " kbp")) +
            scale_x_continuous(breaks=seq(0, 60, by=10)) +
            scale_y_continuous(limits=c(0, ymax),
                               name=expression(pi/D[xy])) +
            scale_colour_manual(values="darkgreen") +
            scale_fill_manual(values="darkgreen") +
            guides(colour=FALSE, fill=FALSE) +
            facet_grid(Species ~ Scaffold,
                       scales="free_x",
                       space="free_x") +
            theme(axis.title=element_text(face="bold"),
                  panel.spacing.x=unit(0, "cm"),
                  text=element_text(size=fontsize),
                  axis.text=element_text(size=fontsize-2),
                  axis.title.y.left=element_text(colour="darkgreen"),
                  axis.ticks.y.left=element_line(colour="darkgreen"),
                  axis.text.y.left=element_text(colour="darkgreen"))
#      print(pi_over_dxy_plot)
      ggsave(paste0('HKA_4fold_', avg_type, 'Avg_w', window_size, '_pi_over_Dxy.pdf'),
             plot=pi_over_dxy_plot,
             width=width_cm,
             height=height_cm,
             units="cm",
             dpi=500)
   }
}
