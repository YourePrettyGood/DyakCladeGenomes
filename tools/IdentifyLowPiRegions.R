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
fai_fn <- options[1]
pi_fn <- options[2]
plot_prefix <- options[3]
chroms <- options[-c(1,2,3)]

#Plotting details:
colorblind_palette <- c("#009E73", "#E79F00", "#9AD0F3", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#000000")

#R functions for finding fuzzy windows:
P_segment <- function(pi, scaffold, threshold) {
   pi %>%
      filter(Scaffold == scaffold) %>%
      mutate(PrefixSum=cumsum(Pi),
             AboveThresh=cumsum(Pi >= threshold)) %>%
      mutate(muP=PrefixSum/row_number(),
             AboveFraction=AboveThresh/row_number())
}
Z_segment <- function(pi, scaffold, threshold) {
   pi %>%
      filter(Scaffold == scaffold) %>%
      arrange(desc(WindowStart)) %>%
      mutate(SuffixSum=cumsum(Pi),
             AboveThresh=cumsum(Pi >= threshold)) %>%
      mutate(muZ=SuffixSum/row_number(),
             AboveFraction=AboveThresh/row_number()) %>%
      arrange(WindowStart)
}
#Single wrapping function to find both P and Z given thresh_[PZ] and max_above_[PZ]:
low_pi_regions <- function(pi, scaffold, thresh_P=NULL, thresh_Z=NULL, max_above_P=1, max_above_Z=1) {
   #The default threshold for pi is an approximation of the average along
   # the selected scaffold:
   #(The approximation is a naive average of window averages.)
   avg_pi <- pi %>%
      mutate(Pi=case_when(is.na(Pi) ~ 0.0,
                          TRUE ~ Pi)) %>%
      filter(Scaffold == scaffold) %>%
      summarize(ApproxPi=sum(Pi)/n());
   if (is.null(thresh_P)) {
      thresh_P <- avg_pi$ApproxPi;
   };
   #The default for the suffix sum threshold is to match the prefix sum threshold:
   if (is.null(thresh_Z)) {
      thresh_Z <- thresh_P;
   };
   #Calculate the prefix cumulative sums:
   #(The mutate mitigates issues caused by NAs)
   Ps <- P_segment(pi %>%
                      mutate(Pi=case_when(is.na(Pi) ~ 0.0,
                                          TRUE ~ Pi)),
                   scaffold, thresh_P);
   #Calculate the suffix cumulative sums:
   #(As above, the mutate mitigates issues caused by NAs)
   Zs <- Z_segment(pi %>%
                      mutate(Pi=case_when(is.na(Pi) ~ 0.0,
                                          TRUE ~ Pi)),
                   scaffold, thresh_Z);
   #Find the best hit for a prefix segment:
   #"Best" here implies the longest segment that has at most max_above_P
   # windows above the pi threshold (thresh_P)
   #We use a simple little trick by checking for over-extension using
   # duplicated() (i.e. the last window in the segment cannot be above
   # the threshold), and find the endpoint by inverting the list and
   # selecting the first element.
   best_P <- Ps %>%
      filter(AboveThresh <= max_above_P) %>%
      filter(duplicated(AboveThresh)) %>%
      arrange(desc(WindowStart)) %>%
      filter(row_number() == 1);
   #Find the best hit for a suffix segment:
   #"Best" here implies the longest segment that has at most max_above_Z
   # windows above the pi threshold (thresh_Z)
   #We use a simple little trick by checking for over-extension using
   # duplicated() (i.e. the first window in the segment cannot be above
   # the threshold), and find the endpoint by uninverting the list and
   # selecting the first element.
   best_Z <- Zs %>%
      filter(AboveThresh <= max_above_Z) %>%
      arrange(desc(WindowStart)) %>%
      filter(duplicated(AboveThresh)) %>%
      arrange(WindowStart) %>%
      filter(row_number() == 1);
   #Output both segments:
   list(P=best_P, Z=best_Z)
}

#Function for generating BED lines from the low_pi_regions() output:
PZ_to_BED <- function(PZ, fai, window_size=100000) {
   P_chrom <- PZ$P$Scaffold
   P_start <- 0
   P_end <- ifelse(PZ$P$WindowStart + window_size - 1 > fai[P_chrom, "Length"],
                   fai[P_chrom, "Length"],
                   PZ$P$WindowStart + window_size - 1)
   Z_chrom <- PZ$Z$Scaffold
   Z_start <- PZ$Z$WindowStart - 1
   Z_end <- fai[Z_chrom, "Length"]
   data.frame(Scaffold=c(P_chrom, Z_chrom),
              BEDStart=c(P_start, Z_start),
              BEDEnd=c(P_end, Z_end))
}

#Function for plotting pi with the low-pi regions highlighted:
low_pi_plot <- function(pi, bed, scaffold, plot_prefix, pi_max=3) {
   #Plotting details:
   linesize <- 0.5
   fontsize <- 12
   keysize <- 30
   ytickfontsize <- 8
   tmargintrim <- 10
   bmargintrim <- 2
   pi_max <- (pi %>%
      filter(Scaffold == scaffold) %>%
      summarize(PiMax=max(Pi, na.rm=TRUE)))$PiMax * 1.05
   #Plot as line with transparent colored rectangles for low-pi,
   # shifting positions into Mbp scaling and scaling pi as percentage:
   pi %>%
      filter(Scaffold == scaffold) %>%
      ggplot(aes(x=WindowStart/1e6, y=Pi*100)) +
         geom_line() +
         geom_rect(data=bed,
                   aes(xmin=(BEDStart+1)/1e6, xmax=BEDEnd/1e6,
                       ymin=0, ymax=pi_max*100),
                   color="red",
                   alpha=0.3,
                   inherit.aes=FALSE) +
         theme_bw() +
         labs(x="Position (Mbp)",
              y=expression(pi~"(%)"),
              title=paste(plot_prefix, chrom, "low pi regions")) +
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
}

#Read in the scaffold lengths from the FAI file:
fai <- read.table(fai_fn,
                  header=FALSE,
                  colClasses=c("character", rep("numeric", 4)),
                  col.names=c("Scaffold", "Length", "Offset", "Wrap", "WrapLF"),
                  sep="\t")
rownames(fai) <- fai$Scaffold

#Read in the windowed pi estimates:
pi <- read.table(gzfile(pi_fn),
                 header=FALSE,
                 colClasses=c("character", rep("numeric", 3)),
                 col.names=c("Scaffold", "WindowStart", "Pi", "UsableFraction"),
                 sep="\t")

#Quick check and fix for the special case of Dsim where the FlyBase chromosome
# arms are annoyingly named with a Scf_ prefix:
#Also adjust for the FlyBase extended scaffold name after whitespace
# (this is an annoying side-effect of my calculatePolymorphism code not
#  trimming scaffold names after the first whitespace)
Scf_add <- ""
if (length(intersect(chroms, fai$Scaffold)) == 0 || length(intersect(chroms, pi$Scaffold)) == 0) {
   if (any(str_detect(fai$Scaffold, "Scf_"))) {
      Scf_add <- "Scf_"
   }
   fai <- fai %>%
      mutate(Scaffold=str_replace(Scaffold, "Scf_", ""))
   rownames(fai) <- fai$Scaffold
   pi <- pi %>%
      mutate(Scaffold=str_replace(Scaffold, "Scf_", "")) %>%
      mutate(Scaffold=str_split_fixed(Scaffold, " ", 2)[,1])
}

#Identify low-pi regions, make QC plots of them, and output a BED of them
# to STDOUT:
pdf(paste0(plot_prefix, "_low_pi.pdf"),
    width=16.0/2.54,
    height=12.0/2.54,
    title=paste0(plot_prefix, " low pi regions"))
for (chrom in chroms) {
   #Check to make sure the chromosome exists in the pi file:
   stopifnot((pi %>% filter(Scaffold == chrom) %>% nrow()) > 0)
   #Determine the prefix and suffix segments for the current chromosome:
   PZ <- low_pi_regions(pi, chrom)
   #Translate the segment lines into BED lines:
   PZ_bed <- PZ_to_BED(PZ, fai)
   #Generate pi plots indicating the low pi regions for QC:
   print(low_pi_plot(pi, PZ_bed, chrom, plot_prefix))
   #Re-add the Scf_ prefix if it was in the input:
   PZ_bed <- PZ_bed %>% mutate(Scaffold=str_c(Scf_add, Scaffold))
   #Output the BED lines:
   writeLines(format_tsv(PZ_bed, col_names=FALSE), stdout())
}
invisible(dev.off())
