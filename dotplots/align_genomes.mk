SHELL=/bin/bash

#Generate pairwise genome alignments with MUMmer
#Assumes refs are available in subdirectories
# of parent directory
#Also assumes that contigs used to compose refs are in the contigs
# subdirectory of the refs directory

#Expects the following in your PATH:
#nucmer
#delta-filter
#mummerplot
#gnuplot
#Additional awk scripts in the ../tools directory:
#mummerPlotOrder.awk
#mummerFixGP.awk
#contigPlotOrder.awk

#Options to adjust:
#Outgroup (prefix of the FASTA for it in the ../refs directory):
OGID := Dmel_ISO1
#Comma-separated (no spaces) list of chromosomes in the order to be plotted:
CHROMS := X,2L,2R,3L,3R,4

#Old reference to align contigs to:
OLDID := Dyak_caf1
#New assemblies to split into contigs and align to old ref:
#Also to align to each other
NEWASMS := Dyak_NY73PB Dyak_Tai18E2

#DO NOT CHANGE THE BELOW:
#Subdirectories that need to be created:
SUBDIRS := deltas logs plot_order plots

#References inferred from directory above the current:
REFS := $(basename $(notdir $(wildcard ../refs/*.fasta)))
#Assemblies that are not the outgroup:
ASMS := $(filter-out $(OGID),$(REFS))
#Assembly-specific FAIs:
FAIS := $(addprefix plot_order/,$(addsuffix .fasta.fai,$(REFS)))
#Assembly-specific plot orders:
REFORDERS := $(addprefix plot_order/,$(addsuffix _scaf_plot_order.tsv,$(REFS) $(OLDID)))
#Contig sets:
CTGS := $(patsubst %_mtDNA_filtered_renamed.fasta,%,$(notdir $(wildcard ../refs/contigs/*.fasta)))
#Contig plot orders:
CTGORDERS := $(addprefix plot_order/,$(addsuffix _contig_plot_order.tsv,$(CTGS)))
#Dotplots of assemblies against the outgroup:
OGPLOTS := $(addprefix plots/$(OGID)_vs_,$(addsuffix _l100.png,$(ASMS)))
#Intermediate GP files:
OGFIXEDGPS := $(addprefix plots/$(OGID)_vs_,$(addsuffix _l100_fixed.gp,$(ASMS)))
OGGPS := $(addprefix plots/$(OGID)_vs_,$(addsuffix _l100.gp,$(ASMS)))
#Intermediate delta files:
OGDELTAS := $(addprefix deltas/$(OGID)_vs_,$(addsuffix _l100.delta,$(ASMS)))
#Dotplots of new assemblies against old:
OLDPLOTS := $(addprefix plots/$(OLDID)_vs_,$(addsuffix _l100c1000.png,$(NEWASMS)))
#Intermediate GP files:
OLDFIXEDGPS := $(addprefix plots/$(OLDID)_vs_,$(addsuffix _l100c1000_fixed.gp,$(NEWASMS)))
OLDGPS := $(addprefix plots/$(OLDID)_vs_,$(addsuffix _l100c1000.gp,$(NEWASMS)))
#Intermediate delta files:
OLDDELTAS := $(addprefix deltas/$(OLDID)_vs_,$(addsuffix _l100c1000.delta,$(NEWASMS)))
#Dotplot between assemblies (e.g. Dyak NY73PB vs. Tai18E2):
#Simplistic variable definition here, only expects length 2 NEWASMS
ASMPLOT := plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.png
#Intermediate GP files:
ASMFIXEDGP := plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000_fixed.gp
ASMGP := plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.gp
#Intermediate delta file:
ASMDELTA := deltas/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.delta
#Dotplots of contigs against their respective chromosomal assemblies:
CTGPLOTS := $(addprefix plots/,$(addsuffix _ctgs_vs_chroms_l100c1000.png,$(CTGS)))
#Intermediate GP files:
CTGFIXEDGPS := $(addprefix plots/,$(addsuffix _ctgs_vs_chroms_l100c1000_fixed.gp,$(CTGS)))
CTGGPS := $(addprefix plots/,$(addsuffix _ctgs_vs_chroms_l100c1000.gp,$(CTGS)))
#Intermediate delta files:
CTGDELTAS := $(addprefix deltas/,$(addsuffix _ctgs_vs_chroms_l100c1000.delta,$(CTGS)))

.PHONY : all vsOG vsOld sameSpecies ctgs clean usage

.SECONDARY : 

usage :
	@echo "Usage:"
	@echo "make -f calculate_pi.mk [task]"
	@echo "Tasks:"
	@echo "all -> Run vsOG, vsOld, sameSpecies, and ctgs"
	@echo "vsOG -> Align new assemblies to outgroup (Dmel)"
	@echo "vsOld -> Align new assembly to old reference (Dyak_caf1)"
	@echo "sameSpecies -> Align new assemblies from the same species"
	@echo "ctgs -> Align contigs against respective chromosomal assemblies"
	@echo "clean -> Clean up all output files and directories"

all : $(OGPLOTS) $(OLDPLOTS) $(ASMPLOT) $(CTGPLOTS)

vsOG : $(OGPLOTS)

vsOld : $(OLDPLOTS)

sameSpecies : $(ASMPLOT)

ctgs : $(CTGPLOTS)

$(OGPLOTS) : plots/$(OGID)_vs_%_l100.png : plots/$(OGID)_vs_%_l100_fixed.gp
	gnuplot $<

$(OLDPLOTS) : plots/$(OLDID)_vs_%_l100c1000.png : plots/$(OLDID)_vs_%_l100c1000_fixed.gp
	gnuplot $<

$(ASMPLOT) : plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000_fixed.gp
	gnuplot $<

$(CTGPLOTS) : plots/%_ctgs_vs_chroms_l100c1000.png : plots/%_ctgs_vs_chroms_l100c1000_fixed.gp
	gnuplot $<

$(OGFIXEDGPS) : plots/$(OGID)_vs_%_l100_fixed.gp : plots/$(OGID)_vs_%_l100.gp
	../tools/mummerFixGP.awk -v "qid=$(OGID)" -v "rid=$*" $< > $@

$(CTGFIXEDGPS) : plots/%_ctgs_vs_chroms_l100c1000_fixed.gp : plots/%_ctgs_vs_chroms_l100c1000.gp
	../tools/mummerFixGP.awk -v "qid=$* contigs" -v "rid=$* arms" $< > $@

$(OLDFIXEDGPS) : plots/$(OLDID)_vs_%_l100c1000_fixed.gp : plots/$(OLDID)_vs_%_l100c1000.gp
	../tools/mummerFixGP.awk -v "qid=$(OLDID)" -v "rid=$*" $< > $@

plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000_fixed.gp : plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.gp
	../tools/mummerFixGP.awk -v "qid=$(word 1,$(NEWASMS))" -v "rid=$(word 2,$(NEWASMS))" $< > $@

$(OGGPS) : plots/$(OGID)_vs_%_l100.gp : plot_order/$(OGID)_scaf_plot_order.tsv plot_order/%_scaf_plot_order.tsv deltas/$(OGID)_vs_%_l100.delta plots
	mummerplot --png --large --filter -Q $(word 1,$^) -R $(word 2,$^) -p plots/$(OGID)_vs_$*_l100 $(word 3,$^) 2> logs/mummerplot_$(OGID)_vs_$*_l100.stderr

$(CTGGPS) : plots/%_ctgs_vs_chroms_l100c1000.gp : plot_order/%_contig_plot_order.tsv plot_order/%_scaf_plot_order.tsv deltas/%_ctgs_vs_chroms_l100c1000.delta plots
	mummerplot --png --large --filter -Q $(word 1,$^) -R $(word 2,$^) -p plots/$*_ctgs_vs_chroms_l100c1000 $(word 3,$^) 2> logs/mummerplot_$*_ctgs_vs_chroms_l100c1000.stderr

$(OLDGPS) : plots/$(OLDID)_vs_%_l100c1000.gp : plot_order/$(OLDID)_scaf_plot_order.tsv plot_order/%_scaf_plot_order.tsv deltas/$(OLDID)_vs_%_l100c1000.delta plots
	mummerplot --png --large --filter -Q $(word 1,$^) -R $(word 2,$^) -p plots/$(OLDID)_vs_$*_l100c1000 $(word 3,$^) 2> logs/mummerplot_$(OLDID)_vs_$*_l100c1000.stderr

plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.gp : plot_order/$(word 1,$(NEWASMS))_scaf_plot_order.tsv plot_order/$(word 2,$(NEWASMS))_scaf_plot_order.tsv deltas/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.delta plots
	mummerplot --png --large --filter -Q $(word 1,$^) -R $(word 2,$^) -p plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000 $(word 3,$^) 2> logs/mummerplot_$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.stderr

$(CTGORDERS) : plot_order/%_contig_plot_order.tsv : ../refs/contigs/%_chromosome_arms_configString.txt plot_order/%_mtDNA_filtered_renamed.fasta.fai
	../tools/contigPlotOrder.awk -v "chromlist=$(CHROMS)" $(word 1,$^) $(word 2,$^) > $@

$(REFORDERS) : plot_order/%_scaf_plot_order.tsv : plot_order/%.fasta.fai
	../tools/mummerPlotOrder.awk -v "chromlist=$(CHROMS)" $< > $@

plot_order/%_mtDNA_filtered_renamed.fasta.fai : ../refs/contigs/%_mtDNA_filtered_renamed.fasta plot_order
	ln -s ../$(word 1,$^) plot_order/$*_mtDNA_filtered_renamed.fasta; samtools faidx plot_order/$*_mtDNA_filtered_renamed.fasta

plot_order/$(OLDID).fasta.fai : ../refs/old/$(OLDID).fasta plot_order
	ln -s ../$(word 1,$^) plot_order/$(OLDID).fasta; samtools faidx plot_order/$(OLDID).fasta

plot_order/%.fasta.fai : ../refs/%.fasta plot_order
	ln -s ../$(word 1,$^) plot_order/$*.fasta; samtools faidx plot_order/$*.fasta

$(OGDELTAS) : deltas/$(OGID)_vs_%_l100.delta : ../refs/%.fasta ../refs/$(OGID).fasta logs deltas
	/usr/bin/time -v nucmer -l100 -p deltas/$(OGID)_vs_$*_l100 $(word 1,$^) $(word 2,$^) 2> logs/nucmer_$(OGID)_vs_$*_l100.stderr > logs/nucmer_$(OGID)_vs_$*_l100.stdout

$(CTGDELTAS) : deltas/%_ctgs_vs_chroms_l100c1000.delta : ../refs/%.fasta ../refs/contigs/%_mtDNA_filtered_renamed.fasta logs deltas
	/usr/bin/time -v nucmer -l100 -c1000 -p deltas/$*_ctgs_vs_chroms_l100c1000 $(word 1,$^) $(word 2,$^) 2> logs/nucmer_$*_ctgs_vs_chroms_l100c1000.stderr > logs/nucmer_$*_ctgs_vs_chroms_l100c1000.stdout

$(OLDDELTAS) : deltas/$(OLDID)_vs_%_l100c1000.delta : ../refs/%.fasta ../refs/old/$(OLDID).fasta logs deltas
	/usr/bin/time -v nucmer -l100 -c1000 -p deltas/$(OLDID)_vs_$*_l100c1000 $(word 1,$^) $(word 2,$^) 2> logs/nucmer_$(OLDID)_vs_$*_l100c1000.stderr > logs/nucmer_$(OLDID)_vs_$*_l100c1000.stdout

deltas/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.delta : ../refs/$(word 1,$(NEWASMS)).fasta ../refs/$(word 2,$(NEWASMS)).fasta logs deltas
	/usr/bin/time -v nucmer -l100 -c1000 -p deltas/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000 $(word 2,$^) $(word 1,$^) 2> logs/nucmer_$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.stderr > logs/nucmer_$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.stdout

$(SUBDIRS) :
	mkdir -p $@

clean :
	for i in $(SUBDIRS); do rm -f $${i}/*.fasta $${i}/*.fai $${i}/*.tsv $${i}/*.delta $${i}/*.[fr]plot $${i}/*.gp $${i}/*.filter $${i}/*.png $${i}/*.stderr $${i}/*.stdout; [[ ! -d $${i} ]] || rmdir $${i}; done
