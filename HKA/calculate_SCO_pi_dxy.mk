SHELL=/bin/bash

#Calculate pi and Dxy for sliding windows based on per-SCO estimates
# from Polymorphorama, including bootstrap estimates to establish
# bootstrap confidence intervals on these estimates.
#The bootstraps can be paired due to usage of the same seeds, thereby
# forming a valid basis for both a per-statistic CI and a ratio CI,
# which we use for the HKA-like statistic: pi/Dxy
#Unfortunately, I couldn't find a decent way to make this makefile
# broadly generalizable for arbitrary refs, primarily due to how we
# specifically require outputs from two different Polymorphrama runs:
# The Dsan, Dtei, and Dyak Polymorphorama runs used Dmel as outgroup
# for Dxy calculations, and the Dmel and Dsim Polymorphorama runs
# used Dtei as outgroup.
#This way we assure equivalence of the divergence metric between the
# two clades, enabling comparison between the clades.
#Assumes refs are available in subdirectories of parent directory,
# and the Polymorphorama output files are available in this directory.

#Expects the following in your PATH:
#bedtools
#samtools
#GNU parallel
#Expects these Awk, C++, and Perl scripts in the tools subdirectory:
#OGwindowMap.awk
#filterWindowOGs.awk
#OGwindowAverage.awk
#bootstrapWindowOGs.awk
#These Awk, C++, and Perl scripts are available in the tools subdirectory
# of the Github repository.
#Also expects the following R script in the tools subdirectory:
#IdentifyLowPiRegions.R

#Options to adjust:
#Window size (bp):
WINDOW := 1000000
#Slide/step (bp) -- preferably an integer multiple of WINDOW:
SLIDE := 200000
#Number of bootstraps:
NBOOTSTRAPS := 10000
#First outgroup:
OUTGROUPONE := Dtei
#Species to evaluate with first outgroup:
CLADEONE := Dmel Dsim
#Second outgroup:
OUTGROUPTWO := Dmel
#Species to evaluate with second outgroup:
CLADETWO := Dsan Dtei Dyak
#Site type to analyze (4fold or FEI):
#I've only actually tested with 4fold
SITETYPE := 4fold
#This variable is only meant for the column names when averaging,
# as Kevin's header uses "4Fold" rather than "4fold":
SITETYPEKD := $(subst f,F,$(SITETYPE))
#References to exclude so that each species is represented by a single
# assembly:
EXCLUDEDREFS := Dyak_Tai18E2
#Polymorphorama compiled and filtered outputs (Kevin prepared these from
# the alignment FASTAs I generated in the CDS and introns directories):
#CLADEONEDATA := Final_filtered_Polymorphorama_dataset_codeml_$(SITETYPE).$(OUTGROUPONE)_root.tsv
CLADEONEDATA := Final_filtered_Polymorphorama_dataset_codeml_$(SITETYPE).$(OUTGROUPONE)_root.low_pi_included.tsv
#CLADETWODATA := Final_filtered_Polymorphorama_dataset_codeml_$(SITETYPE).$(OUTGROUPTWO)_root.tsv
CLADETWODATA := Final_filtered_Polymorphorama_dataset_codeml_$(SITETYPE).$(OUTGROUPTWO)_root.low_pi_included.tsv
#Header names for the orthogroup and species columns in the above files:
OGCOL := ortho
SPPCOL := pop

#DO NOT CHANGE THE BELOW:
#References inferred from directory above the current:
REFS := $(filter-out $(EXCLUDEDREFS),$(basename $(notdir $(wildcard ../refs/*.fasta))))
#Subdirectories that need to be created:
SUBDIRS := BEDs maps stats/bootstraps stats/observed logs

#Rdata file for plotting and any other analyses:
RDATA := HKA_$(SITETYPE)_stats_wCIs.Rdata
#Offsets file for lining up windows between species:
OFFSETS := Scaffold_OG_offsets.tsv
#Window BED files:
WINDOWBEDS := $(addprefix BEDs/,$(addsuffix _w$(WINDOW)_s$(SLIDE)_windows.bed,$(CLADEONE) $(CLADETWO)))
#Orthogroup-to-window maps:
WINDOWMAPS := $(addprefix maps/,$(addsuffix _w$(WINDOW)_s$(SLIDE)_windowmap.tcsv,$(CLADEONE) $(CLADETWO)))
#Observed average pi and Dxy:
CLADEONENOBSDXY := $(addsuffix _Dxy_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/observed/,$(foreach spp,$(CLADEONE),$(spp)_w$(WINDOW)_s$(SLIDE))))
CLADEONEWOBSDXY := $(addsuffix _Dxy_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/observed/,$(foreach spp,$(CLADEONE),$(spp)_w$(WINDOW)_s$(SLIDE))))
CLADEONENOBSPI := $(addsuffix _Pi_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/observed/,$(foreach spp,$(CLADEONE),$(spp)_w$(WINDOW)_s$(SLIDE))))
CLADEONEWOBSPI := $(addsuffix _Pi_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/observed/,$(foreach spp,$(CLADEONE),$(spp)_w$(WINDOW)_s$(SLIDE))))
CLADETWONOBSDXY := $(addsuffix _Dxy_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/observed/,$(foreach spp,$(CLADETWO),$(spp)_w$(WINDOW)_s$(SLIDE))))
CLADETWOWOBSDXY := $(addsuffix _Dxy_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/observed/,$(foreach spp,$(CLADETWO),$(spp)_w$(WINDOW)_s$(SLIDE))))
CLADETWONOBSPI := $(addsuffix _Pi_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/observed/,$(foreach spp,$(CLADETWO),$(spp)_w$(WINDOW)_s$(SLIDE))))
CLADETWOWOBSPI := $(addsuffix _Pi_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/observed/,$(foreach spp,$(CLADETWO),$(spp)_w$(WINDOW)_s$(SLIDE))))
#Window to orthogroup maps intended for bootstrapping:
CLADEONEOGDXYMAPS := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv,$(addprefix maps/,$(CLADEONE)))
CLADEONEOGPIMAPS := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv,$(addprefix maps/,$(CLADEONE)))
CLADETWOOGDXYMAPS := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv,$(addprefix maps/,$(CLADETWO)))
CLADETWOOGPIMAPS := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv,$(addprefix maps/,$(CLADETWO)))
#Bootstrapped average pi and Dxy:
#I had to hard-code this because I couldn't figure out how to do a two-part
# pattern rule (e.g. a single pattern rule that works for all the Dxy
# bootstraps and works with parallelization). For now, the pattern rule
# abstracts out the bootstrap seed, but we hard-code the species, statistic,
# and averaging method. In an ideal world, we would abstract across all
# four of these variables with a single pattern rule for the bootstraps.
#If anyone reading this knows of a way to achieve this in GNU make,
# please let me know! (I realize stuff like snakemake and nextflow could
# do this.)
DMELNBOOTDXY := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/bootstraps/Dmel_,$(shell seq 1 $(NBOOTSTRAPS))))
DMELWBOOTDXY := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/bootstraps/Dmel_,$(shell seq 1 $(NBOOTSTRAPS))))
DMELNBOOTPI := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/bootstraps/Dmel_,$(shell seq 1 $(NBOOTSTRAPS))))
DMELWBOOTPI := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/bootstraps/Dmel_,$(shell seq 1 $(NBOOTSTRAPS))))
DSIMNBOOTDXY := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/bootstraps/Dsim_,$(shell seq 1 $(NBOOTSTRAPS))))
DSIMWBOOTDXY := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/bootstraps/Dsim_,$(shell seq 1 $(NBOOTSTRAPS))))
DSIMNBOOTPI := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/bootstraps/Dsim_,$(shell seq 1 $(NBOOTSTRAPS))))
DSIMWBOOTPI := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/bootstraps/Dsim_,$(shell seq 1 $(NBOOTSTRAPS))))
DSANNBOOTDXY := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/bootstraps/Dsan_,$(shell seq 1 $(NBOOTSTRAPS))))
DSANWBOOTDXY := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/bootstraps/Dsan_,$(shell seq 1 $(NBOOTSTRAPS))))
DSANNBOOTPI := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/bootstraps/Dsan_,$(shell seq 1 $(NBOOTSTRAPS))))
DSANWBOOTPI := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/bootstraps/Dsan_,$(shell seq 1 $(NBOOTSTRAPS))))
DTEINBOOTDXY := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/bootstraps/Dtei_,$(shell seq 1 $(NBOOTSTRAPS))))
DTEIWBOOTDXY := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/bootstraps/Dtei_,$(shell seq 1 $(NBOOTSTRAPS))))
DTEINBOOTPI := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/bootstraps/Dtei_,$(shell seq 1 $(NBOOTSTRAPS))))
DTEIWBOOTPI := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/bootstraps/Dtei_,$(shell seq 1 $(NBOOTSTRAPS))))
DYAKNBOOTDXY := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/bootstraps/Dyak_,$(shell seq 1 $(NBOOTSTRAPS))))
DYAKWBOOTDXY := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/bootstraps/Dyak_,$(shell seq 1 $(NBOOTSTRAPS))))
DYAKNBOOTPI := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz,$(addprefix stats/bootstraps/Dyak_,$(shell seq 1 $(NBOOTSTRAPS))))
DYAKWBOOTPI := $(addsuffix _w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz,$(addprefix stats/bootstraps/Dyak_,$(shell seq 1 $(NBOOTSTRAPS))))
#Combined files for both observed and bootstrap averages:
DMELNALLDXY := stats/Dmel_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg_obs_wBootstraps.tsv.gz
DMELWALLDXY := stats/Dmel_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg_obs_wBootstraps.tsv.gz
DMELNALLPI := stats/Dmel_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg_obs_wBootstraps.tsv.gz
DMELWALLPI := stats/Dmel_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg_obs_wBootstraps.tsv.gz
DSIMNALLDXY := stats/Dsim_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg_obs_wBootstraps.tsv.gz
DSIMWALLDXY := stats/Dsim_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg_obs_wBootstraps.tsv.gz
DSIMNALLPI := stats/Dsim_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg_obs_wBootstraps.tsv.gz
DSIMWALLPI := stats/Dsim_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg_obs_wBootstraps.tsv.gz
DSANNALLDXY := stats/Dsan_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg_obs_wBootstraps.tsv.gz
DSANWALLDXY := stats/Dsan_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg_obs_wBootstraps.tsv.gz
DSANNALLPI := stats/Dsan_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg_obs_wBootstraps.tsv.gz
DSANWALLPI := stats/Dsan_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg_obs_wBootstraps.tsv.gz
DTEINALLDXY := stats/Dtei_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg_obs_wBootstraps.tsv.gz
DTEIWALLDXY := stats/Dtei_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg_obs_wBootstraps.tsv.gz
DTEINALLPI := stats/Dtei_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg_obs_wBootstraps.tsv.gz
DTEIWALLPI := stats/Dtei_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg_obs_wBootstraps.tsv.gz
DYAKNALLDXY := stats/Dyak_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg_obs_wBootstraps.tsv.gz
DYAKWALLDXY := stats/Dyak_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg_obs_wBootstraps.tsv.gz
DYAKNALLPI := stats/Dyak_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg_obs_wBootstraps.tsv.gz
DYAKWALLPI := stats/Dyak_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg_obs_wBootstraps.tsv.gz
#A single variable to act as the target list for each species:
DMELSTATS := $(DMELNALLDXY) $(DMELWALLDXY) $(DMELNALLPI) $(DMELWALLPI)
DSIMSTATS := $(DSIMNALLDXY) $(DSIMWALLDXY) $(DSIMNALLPI) $(DSIMWALLPI)
DSANSTATS := $(DSANNALLDXY) $(DSANWALLDXY) $(DSANNALLPI) $(DSANWALLPI)
DTEISTATS := $(DTEINALLDXY) $(DTEIWALLDXY) $(DTEINALLPI) $(DTEIWALLPI)
DYAKSTATS := $(DYAKNALLDXY) $(DYAKWALLDXY) $(DYAKNALLPI) $(DYAKWALLPI)

.PHONY : all clean usage

.SECONDARY : 

usage :
	@echo "Usage:"
	@echo "make -f calculate_pi.mk [task]"
	@echo "Tasks:"
	@echo "all -> Calculate pi and Dxy with bootstraps from $(SITETYPE) sites in SCOs"
	@echo "clean -> Clean up all output files and directories"

all : $(RDATA)

#Run the R script to collate all the windowed statistics and bootstraps,
# and calculate the 95% CI bounds for each window:
$(RDATA) : $(OFFSETS) $(DMELSTATS) $(DSIMSTATS) $(DSANSTATS) $(DTEISTATS) $(DYAKSTATS)
	@echo "Calculating the 95% CI bounds and plotting pi and Dxy"
	./DyakClade_HKA_4fold.R $(word 1,$^) $(SLIDE) $(WINDOW)

#Determine the offsets required to line up the major arms:
#This is based on the SCO nearest the telomere of each arm
$(OFFSETS) : ../genomewide_pi/orthogroup_locations
	@echo "Determining offsets for lining up major arms"
	sort -k1,1 -k2,2n -k3,3n < $< | ../tools/lineupArms.awk > $@

#Combine the observed and bootstrap files together:
$(DMELNALLDXY) : stats/observed/Dmel_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz $(DMELNBOOTDXY)
	@echo "Compiling observed Dxy (naive average) and bootstraps for Dmel"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dmel_$${i}_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz; done | gzip -9 >> $@

$(DMELWALLDXY) : stats/observed/Dmel_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz $(DMELWBOOTDXY)
	@echo "Compiling observed Dxy (weighted average) and bootstraps for Dmel"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dmel_$${i}_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz; done | gzip -9 >> $@

$(DMELNALLPI) : stats/observed/Dmel_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz $(DMELNBOOTPI)
	@echo "Compiling observed Pi (naive average) and bootstraps for Dmel"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dmel_$${i}_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz; done | gzip -9 >> $@

$(DMELWALLPI) : stats/observed/Dmel_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz $(DMELWBOOTPI)
	@echo "Compiling observed Pi (weighted average) and bootstraps for Dmel"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dmel_$${i}_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz; done | gzip -9 >> $@

$(DSIMNALLDXY) : stats/observed/Dsim_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz $(DSIMNBOOTDXY)
	@echo "Compiling observed Dxy (naive average) and bootstraps for Dsim"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dsim_$${i}_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz; done | gzip -9 >> $@

$(DSIMWALLDXY) : stats/observed/Dsim_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz $(DSIMWBOOTDXY)
	@echo "Compiling observed Dxy (weighted average) and bootstraps for Dsim"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dsim_$${i}_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz; done | gzip -9 >> $@

$(DSIMNALLPI) : stats/observed/Dsim_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz $(DSIMNBOOTPI)
	@echo "Compiling observed Pi (naive average) and bootstraps for Dsim"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dsim_$${i}_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz; done | gzip -9 >> $@

$(DSIMWALLPI) : stats/observed/Dsim_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz $(DSIMWBOOTPI)
	@echo "Compiling observed Pi (weighted average) and bootstraps for Dsim"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dsim_$${i}_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz; done | gzip -9 >> $@

$(DSANNALLDXY) : stats/observed/Dsan_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz $(DSANNBOOTDXY)
	@echo "Compiling observed Dxy (naive average) and bootstraps for Dsan"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dsan_$${i}_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz; done | gzip -9 >> $@

$(DSANWALLDXY) : stats/observed/Dsan_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz $(DSANWBOOTDXY)
	@echo "Compiling observed Dxy (weighted average) and bootstraps for Dsan"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dsan_$${i}_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz; done | gzip -9 >> $@

$(DSANNALLPI) : stats/observed/Dsan_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz $(DSANNBOOTPI)
	@echo "Compiling observed Pi (naive average) and bootstraps for Dsan"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dsan_$${i}_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz; done | gzip -9 >> $@

$(DSANWALLPI) : stats/observed/Dsan_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz $(DSANWBOOTPI)
	@echo "Compiling observed Pi (weighted average) and bootstraps for Dsan"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dsan_$${i}_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz; done | gzip -9 >> $@

$(DTEINALLDXY) : stats/observed/Dtei_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz $(DTEINBOOTDXY)
	@echo "Compiling observed Dxy (naive average) and bootstraps for Dtei"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dtei_$${i}_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz; done | gzip -9 >> $@

$(DTEIWALLDXY) : stats/observed/Dtei_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz $(DTEIWBOOTDXY)
	@echo "Compiling observed Dxy (weighted average) and bootstraps for Dtei"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dtei_$${i}_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz; done | gzip -9 >> $@

$(DTEINALLPI) : stats/observed/Dtei_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz $(DTEINBOOTPI)
	@echo "Compiling observed Pi (naive average) and bootstraps for Dtei"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dtei_$${i}_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz; done | gzip -9 >> $@

$(DTEIWALLPI) : stats/observed/Dtei_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz $(DTEIWBOOTPI)
	@echo "Compiling observed Pi (weighted average) and bootstraps for Dtei"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dtei_$${i}_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz; done | gzip -9 >> $@

$(DYAKNALLDXY) : stats/observed/Dyak_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz $(DYAKNBOOTDXY)
	@echo "Compiling observed Dxy (naive average) and bootstraps for Dyak"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dyak_$${i}_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz; done | gzip -9 >> $@

$(DYAKWALLDXY) : stats/observed/Dyak_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz $(DYAKWBOOTDXY)
	@echo "Compiling observed Dxy (weighted average) and bootstraps for Dyak"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dyak_$${i}_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz; done | gzip -9 >> $@

$(DYAKNALLPI) : stats/observed/Dyak_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz $(DYAKNBOOTPI)
	@echo "Compiling observed Pi (naive average) and bootstraps for Dyak"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dyak_$${i}_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz; done | gzip -9 >> $@

$(DYAKWALLPI) : stats/observed/Dyak_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz $(DYAKWBOOTPI)
	@echo "Compiling observed Pi (weighted average) and bootstraps for Dyak"
	cat $(word 1,$^) > $@
	for i in {1..$(NBOOTSTRAPS)}; do gzip -dc stats/bootstraps/Dyak_$${i}_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz; done | gzip -9 >> $@

#Do the bootstraps:
#We have to hard-code species, statistic, and average type for these due to
# limitations on pattern rules.
#Frustratingly, the column names aren't particularly consistently formatted,
# so 4fold is termed 4Fold, and Pi is pi, but Dxy is Dxy...
$(DMELNBOOTDXY) : stats/bootstraps/Dmel_%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz : maps/Dmel_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv $(CLADEONEDATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dmel" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DMELWBOOTDXY) : stats/bootstraps/Dmel_%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz : maps/Dmel_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv $(CLADEONEDATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dmel" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DMELNBOOTPI) : stats/bootstraps/Dmel_%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz : maps/Dmel_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv $(CLADEONEDATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dmel" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DMELWBOOTPI) : stats/bootstraps/Dmel_%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz : maps/Dmel_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv $(CLADEONEDATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dmel" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DSIMNBOOTDXY) : stats/bootstraps/Dsim_%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz : maps/Dsim_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv $(CLADEONEDATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dsim" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DSIMWBOOTDXY) : stats/bootstraps/Dsim_%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz : maps/Dsim_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv $(CLADEONEDATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dsim" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DSIMNBOOTPI) : stats/bootstraps/Dsim_%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz : maps/Dsim_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv $(CLADEONEDATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dsim" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DSIMWBOOTPI) : stats/bootstraps/Dsim_%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz : maps/Dsim_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv $(CLADEONEDATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dsim" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DSANNBOOTDXY) : stats/bootstraps/Dsan_%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz : maps/Dsan_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv $(CLADETWODATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dsan" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DSANWBOOTDXY) : stats/bootstraps/Dsan_%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz : maps/Dsan_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv $(CLADETWODATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dsan" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DSANNBOOTPI) : stats/bootstraps/Dsan_%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz : maps/Dsan_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv $(CLADETWODATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dsan" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DSANWBOOTPI) : stats/bootstraps/Dsan_%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz : maps/Dsan_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv $(CLADETWODATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dsan" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DTEINBOOTDXY) : stats/bootstraps/Dtei_%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz : maps/Dtei_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv $(CLADETWODATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dtei" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DTEIWBOOTDXY) : stats/bootstraps/Dtei_%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz : maps/Dtei_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv $(CLADETWODATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dtei" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DTEINBOOTPI) : stats/bootstraps/Dtei_%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz : maps/Dtei_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv $(CLADETWODATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dtei" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DTEIWBOOTPI) : stats/bootstraps/Dtei_%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz : maps/Dtei_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv $(CLADETWODATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dtei" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DYAKNBOOTDXY) : stats/bootstraps/Dyak_%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz : maps/Dyak_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv $(CLADETWODATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dyak" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DYAKWBOOTDXY) : stats/bootstraps/Dyak_%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz : maps/Dyak_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv $(CLADETWODATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dyak" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DYAKNBOOTPI) : stats/bootstraps/Dyak_%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz : maps/Dyak_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv $(CLADETWODATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dyak" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

$(DYAKWBOOTPI) : stats/bootstraps/Dyak_%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz : maps/Dyak_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv $(CLADETWODATA) stats/bootstraps
	../tools/bootstrapWindowOGs.awk -v "seed=$*" $(word 1,$^) | ../tools/OGwindowAverage.awk -v "species=Dyak" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=$*" -v "noheader=1" - $(word 2,$^) | gzip -9 > $@

#Get the observed averages:
#Frustratingly, the column names aren't particularly consistently formatted,
# so 4fold is termed 4Fold, and Pi is pi, but Dxy is Dxy...
$(CLADEONENOBSDXY) : stats/observed/%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz : maps/%_w$(WINDOW)_s$(SLIDE)_windowmap.tcsv $(CLADEONEDATA) stats/observed
	@echo "Taking naive average Dxy for each window for $*"
	../tools/OGwindowAverage.awk -v "species=$*" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "id=0" $(word 1,$^) $(word 2,$^) | gzip -9 > $@

$(CLADEONEWOBSDXY) : stats/observed/%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz : maps/%_w$(WINDOW)_s$(SLIDE)_windowmap.tcsv $(CLADEONEDATA) stats/observed
	@echo "Taking weighted average Dxy for each window for $* (weight is number of sites in alignment)"
	../tools/OGwindowAverage.awk -v "species=$*" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=0" $(word 1,$^) $(word 2,$^) | gzip -9 > $@

$(CLADEONENOBSPI) : stats/observed/%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz : maps/%_w$(WINDOW)_s$(SLIDE)_windowmap.tcsv $(CLADEONEDATA) stats/observed
	@echo "Taking naive average Pi for each window for $*"
	../tools/OGwindowAverage.awk -v "species=$*" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "id=0" $(word 1,$^) $(word 2,$^) | gzip -9 > $@

$(CLADEONEWOBSPI) : stats/observed/%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz : maps/%_w$(WINDOW)_s$(SLIDE)_windowmap.tcsv $(CLADEONEDATA) stats/observed
	@echo "Taking weighted average Pi for each window for $* (weight is number of sites in alignment)"
	../tools/OGwindowAverage.awk -v "species=$*" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=0" $(word 1,$^) $(word 2,$^) | gzip -9 > $@

$(CLADETWONOBSDXY) : stats/observed/%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_naiveAvg.tsv.gz : maps/%_w$(WINDOW)_s$(SLIDE)_windowmap.tcsv $(CLADETWODATA) stats/observed
	@echo "Taking naive average Dxy for each window for $*"
	../tools/OGwindowAverage.awk -v "species=$*" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "id=0" $(word 1,$^) $(word 2,$^) | gzip -9 > $@

$(CLADETWOWOBSDXY) : stats/observed/%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_weightedAvg.tsv.gz : maps/%_w$(WINDOW)_s$(SLIDE)_windowmap.tcsv $(CLADETWODATA) stats/observed
	@echo "Taking weighted average Dxy for each window for $* (weight is number of sites in alignment)"
	../tools/OGwindowAverage.awk -v "species=$*" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=0" $(word 1,$^) $(word 2,$^) | gzip -9 > $@

$(CLADETWONOBSPI) : stats/observed/%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_naiveAvg.tsv.gz : maps/%_w$(WINDOW)_s$(SLIDE)_windowmap.tcsv $(CLADETWODATA) stats/observed
	@echo "Taking naive average Pi for each window for $*"
	../tools/OGwindowAverage.awk -v "species=$*" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "id=0" $(word 1,$^) $(word 2,$^) | gzip -9 > $@

$(CLADETWOWOBSPI) : stats/observed/%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_weightedAvg.tsv.gz : maps/%_w$(WINDOW)_s$(SLIDE)_windowmap.tcsv $(CLADETWODATA) stats/observed
	@echo "Taking weighted average Pi for each window for $* (weight is number of sites in alignment)"
	../tools/OGwindowAverage.awk -v "species=$*" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=pi_$(SITETYPEKD)" -v "weightcol=X$(SITETYPEKD)Sites" -v "id=0" $(word 1,$^) $(word 2,$^) | gzip -9 > $@

#Generate the map that maps windows to lists of single-copy orthogroups:
#This map facilitates bootstrapping within windows, since we just split
# the list on a given line (i.e. a given window) and bootstrap that list.
#This won't be a perfect inverse map of the orthogroup-to-window map, as
# we do some filtering out of orthogroups with NA in the statistic of
# interest, since these are useless for bootstrapping.
$(CLADEONEOGDXYMAPS) : maps/%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv : maps/%_w$(WINDOW)_s$(SLIDE)_windowmap.tcsv $(CLADEONEDATA) maps
	@echo "Generating the window-to-orthogroup map for bootstrapping for Dxy of $*"
	../tools/filterWindowOGs.awk -v "species=$*" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" $(word 1,$^) $(word 2,$^) > $@

$(CLADEONEOGPIMAPS) : maps/%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv : maps/%_w$(WINDOW)_s$(SLIDE)_windowmap.tcsv $(CLADEONEDATA) maps
	@echo "Generating the window-to-orthogroup map for bootstrapping for Pi of $*"
	../tools/filterWindowOGs.awk -v "species=$*" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Pi_$(SITETYPEKD)" $(word 1,$^) $(word 2,$^) > $@

$(CLADETWOOGDXYMAPS) : maps/%_w$(WINDOW)_s$(SLIDE)_Dxy_$(SITETYPE)_ogmap.tcsv : maps/%_w$(WINDOW)_s$(SLIDE)_windowmap.tcsv $(CLADETWODATA) maps
	@echo "Generating the window-to-orthogroup map for bootstrapping for Dxy of $*"
	../tools/filterWindowOGs.awk -v "species=$*" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Dxy_$(SITETYPEKD)" $(word 1,$^) $(word 2,$^) > $@

$(CLADETWOOGPIMAPS) : maps/%_w$(WINDOW)_s$(SLIDE)_Pi_$(SITETYPE)_ogmap.tcsv : maps/%_w$(WINDOW)_s$(SLIDE)_windowmap.tcsv $(CLADETWODATA) maps
	@echo "Generating the window-to-orthogroup map for bootstrapping for Pi of $*"
	../tools/filterWindowOGs.awk -v "species=$*" -v "ogcol=$(OGCOL)" -v "sppcol=$(SPPCOL)" -v "statcol=Pi_$(SITETYPEKD)" $(word 1,$^) $(word 2,$^) > $@

#Generate the map that maps single-copy orthogroups to windows:
#This map allows us to easily take the average within each window given
# a table of statistics, one row per orthogroup.
$(WINDOWMAPS) : maps/%_w$(WINDOW)_s$(SLIDE)_windowmap.tcsv : BEDs/%_w$(WINDOW)_s$(SLIDE)_windows.bed ../genomewide_pi/orthogroup_locations maps
	@echo "Generating the map of single-copy orthogroups to windows for $*"
	../tools/OGwindowMap.awk -v "species=$*" $(word 1,$^) $(word 2,$^) > $@

#Generate windows for each species:
$(WINDOWBEDS) : BEDs/%_w$(WINDOW)_s$(SLIDE)_windows.bed : FAIs/%.genome BEDs
	@echo "Generating $(WINDOW) bp windows with $(SLIDE) bp step for $*"
	bedtools makewindows -g $(word 1,$^) -w $(WINDOW) -s $(SLIDE) > $@

FAIs/%.genome : FAIs/%.fasta.fai
	cut -f1,2 $< > $@

FAIs/%.fasta.fai : FAIs/%.fasta
	@echo "Indexing reference $* with samtools"
	samtools faidx $<

FAIs/%.fasta : FAIs
	cd FAIs/; \
	ln -sf $*_*.fasta $*.fasta; \
	cd ..

FAIs : 
	mkdir -p $@; \
	ln -sf $(addprefix ../../refs/,$(addsuffix .fasta,$(REFS))) FAIs/

$(SUBDIRS) :
	mkdir -p $@

clean :
	rm -rf stats/
	for i in $(SUBDIRS); do rm -f $${i}/*.fasta $${i}/*.fai $${i}/*.genome $${i}/*.bed $${i}/*.tcsv $${i}/*.stderr $${i}/*.tsv.gz; [[ ! -d $${i} ]] || rmdir $${i}; done
	rm -rf FAIs/
	rm -f $(OFFSETS) $(RDATA) *.pdf
