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
#New assemblies to align to old ref:
#Also to align to each other and to closely-related species (i.e. Dsan)
NEWASMS := Dyak_NY73PB Dyak_Tai18E2

#Closely-related species to NEWASMS (i.e. Dsan):
CLOSEID := Dsan_STOCAGO1482

#Miller long-read assembly:
MILLER := Dyak_Miller
#New assembly to compare to Miller assembly:
MILLERMATCH := Dyak_Tai18E2

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
REFORDERS := $(addprefix plot_order/,$(addsuffix _scaf_plot_order.tsv,$(REFS) $(OLDID) $(MILLER)))
#Contig sets:
CTGS := $(patsubst %_contigs_filtered_softmasked_w60.fasta,%,$(notdir $(wildcard ../refs/contigs/*_contigs_filtered_softmasked_w60.fasta)))
#Contig plot orders:
CTGORDERS := $(addprefix plot_order/,$(addsuffix _contig_plot_order.tsv,$(CTGS)))
#Unscaffolded contig plot orders:
UNSCAFORDERS := $(addprefix plot_order/,$(addsuffix _unscaf_plot_order.tsv,$(CTGS)))

#Dotplots of scaffolded assemblies against the outgroup:
OGPLOTS := $(addprefix plots/$(OGID)_vs_,$(addsuffix _l100.png,$(ASMS)))
#Intermediate GP files for scaffolds vs. OG:
OGFIXEDGPS := $(addprefix plots/$(OGID)_vs_,$(addsuffix _l100_fixed.gp,$(ASMS)))
OGGPS := $(addprefix plots/$(OGID)_vs_,$(addsuffix _l100.gp,$(ASMS)))
#Intermediate delta files for scaffolds vs. OG:
OGDELTAS := $(addprefix deltas/$(OGID)_vs_,$(addsuffix _l100.delta,$(ASMS)))
#Dotplots of contigs from assemblies against the outgroup:
OGCTGPLOTS := $(addprefix plots/$(OGID)_vs_,$(addsuffix _ctgs_l100.png,$(CTGS)))
#Intermediate GP files for contigs vs. OG:
OGCTGFIXEDGPS := $(addprefix plots/$(OGID)_vs_,$(addsuffix _ctgs_l100_fixed.gp,$(CTGS)))
OGCTGGPS := $(addprefix plots/$(OGID)_vs_,$(addsuffix _ctgs_l100.gp,$(CTGS)))
#Intermediate delta files for contigs vs. OG:
OGCTGDELTAS := $(addprefix deltas/$(OGID)_vs_,$(addsuffix _ctgs_l100.delta,$(CTGS)))

#Dotplots of new assemblies against old:
OLDPLOTS := $(addprefix plots/$(OLDID)_vs_,$(addsuffix _l100c1000.png,$(NEWASMS)))
#Intermediate GP files:
OLDFIXEDGPS := $(addprefix plots/$(OLDID)_vs_,$(addsuffix _l100c1000_fixed.gp,$(NEWASMS)))
OLDGPS := $(addprefix plots/$(OLDID)_vs_,$(addsuffix _l100c1000.gp,$(NEWASMS)))
#Intermediate delta files:
OLDDELTAS := $(addprefix deltas/$(OLDID)_vs_,$(addsuffix _l100c1000.delta,$(NEWASMS)))

#Dotplots of Miller assembly against old and Tai18E2:
MILLERPLOTS := $(addprefix plots/$(MILLER)_vs_,$(addsuffix _l100c1000.png,$(OLDID) $(MILLERMATCH)))
#Intermediate GP files:
MILLERFIXEDGPS := $(addprefix plots/$(MILLER)_vs_,$(addsuffix _l100c1000_fixed.gp,$(OLDID) $(MILLERMATCH)))
MILLERGPS := $(addprefix plots/$(MILLER)_vs_,$(addsuffix _l100c1000.gp,$(OLDID) $(MILLERMATCH)))
#Intermediate delta files:
MILLERDELTAS := $(addprefix deltas/$(MILLER)_vs_,$(addsuffix _l100c1000.delta,$(OLDID) $(MILLERMATCH)))

#Dotplot between assemblies (e.g. Dyak NY73PB vs. Tai18E2):
#Simplistic variable definition here, only expects length 2 NEWASMS
ASMPLOT := plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.png
#Intermediate GP files:
ASMFIXEDGP := plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000_fixed.gp
ASMGP := plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.gp
#Intermediate delta file:
ASMDELTA := deltas/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.delta

#Dotplots between Dyak assemblies and closely-related species (i.e. Dsan):
CLOSEPLOTS := $(addprefix plots/$(CLOSEID)_vs_,$(addsuffix _l100c1000.png,$(NEWASMS)))
#Intermediate GP files:
CLOSEFIXEDGPS := $(addprefix plots/$(CLOSEID)_vs_,$(addsuffix _l100c1000_fixed.gp,$(NEWASMS)))
CLOSEGPS := $(addprefix plots/$(CLOSEID)_vs_,$(addsuffix _l100c1000.gp,$(NEWASMS)))
#Intermediate delta files:
CLOSEDELTAS := $(addprefix deltas/$(CLOSEID)_vs_,$(addsuffix _l100c1000.delta,$(NEWASMS)))

#Dotplots of scaffolded contigs against their respective chromosomal assemblies:
CTGPLOTS := $(addprefix plots/,$(addsuffix _ctgs_vs_chroms_l100c1000.png,$(CTGS)))
#Intermediate GP files:
CTGFIXEDGPS := $(addprefix plots/,$(addsuffix _ctgs_vs_chroms_l100c1000_fixed.gp,$(CTGS)))
CTGGPS := $(addprefix plots/,$(addsuffix _ctgs_vs_chroms_l100c1000.gp,$(CTGS)))
#Intermediate delta files:
CTGDELTAS := $(addprefix deltas/,$(addsuffix _ctgs_vs_chroms_l100c1000.delta,$(CTGS)))

#Dotplots of unscaffolded contigs against their respective chromosomal assemblies:
UNSCAFPLOTS := $(addprefix plots/,$(addsuffix _unscaf_vs_chroms_l100c1000.png,$(CTGS)))
#Intermediate GP files:
UNSCAFFIXEDGPS := $(addprefix plots/,$(addsuffix _unscaf_vs_chroms_l100c1000_fixed.gp,$(CTGS)))
UNSCAFGPS := $(addprefix plots/,$(addsuffix _unscaf_vs_chroms_l100c1000.gp,$(CTGS)))

.PHONY : all vsOG ctgsvsOG vsOld vsMiller sameSpecies vsClose ctgs unscaf clean usage

.SECONDARY : 

usage :
	@echo "Usage:"
	@echo "make -f align_genomes.mk [task]"
	@echo "Tasks:"
	@echo "all -> Run vsOG, vsOld, sameSpecies, ctgs, and unscaf"
	@echo "vsOG -> Align new assemblies to outgroup (Dmel)"
	@echo "ctgsvsOG -> Align contigs of new assemblies to outgroup (Dmel)"
	@echo "vsOld -> Align new assembly to old reference (Dyak_caf1)"
	@echo "vsMiller -> Align Miller assembly to old and new references"
	@echo "            (Dyak_Miller)"
	@echo "sameSpecies -> Align new assemblies from the same species (Dyak)"
	@echo "vsClose -> Align new assemblies to closely-related species (Dsan)"
	@echo "ctgs -> Align scaffolded contigs against respective chromosomal assemblies"
	@echo "        (This provides information about joins)"
	@echo "unscaf -> Align unscaffolded contigs against respective chromosomal assemblies"
	@echo "          (This checks for alternate haplotype contigs)"
	@echo "clean -> Clean up all output files and directories"

all : $(OGPLOTS) $(OGCTGPLOTS) $(OLDPLOTS) $(MILLERPLOTS) $(ASMPLOT) $(CLOSEPLOTS) $(CTGPLOTS) $(UNSCAFPLOTS)

vsOG : $(OGPLOTS)

ctgsvsOG : $(OGCTGPLOTS)

vsOld : $(OLDPLOTS)

vsMiller : $(MILLERPLOTS)

sameSpecies : $(ASMPLOT)

vsClose : $(CLOSEPLOTS)

ctgs : $(CTGPLOTS)

unscaf : $(UNSCAFPLOTS)

$(OGCTGPLOTS) : plots/$(OGID)_vs_%_ctgs_l100.png : plots/$(OGID)_vs_%_ctgs_l100_fixed.gp
	gnuplot $<

$(OGPLOTS) : plots/$(OGID)_vs_%_l100.png : plots/$(OGID)_vs_%_l100_fixed.gp
	gnuplot $<

$(OLDPLOTS) : plots/$(OLDID)_vs_%_l100c1000.png : plots/$(OLDID)_vs_%_l100c1000_fixed.gp
	gnuplot $<

$(MILLERPLOTS) : plots/$(MILLER)_vs_%_l100c1000.png : plots/$(MILLER)_vs_%_l100c1000_fixed.gp
	gnuplot $<

$(ASMPLOT) : plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000_fixed.gp
	gnuplot $<

$(CLOSEPLOTS) : plots/$(CLOSEID)_vs_%_l100c1000.png : plots/$(CLOSEID)_vs_%_l100c1000_fixed.gp
	gnuplot $<

$(CTGPLOTS) : plots/%_ctgs_vs_chroms_l100c1000.png : plots/%_ctgs_vs_chroms_l100c1000_fixed.gp
	gnuplot $<

$(UNSCAFPLOTS) : plots/%_unscaf_vs_chroms_l100c1000.png : plots/%_unscaf_vs_chroms_l100c1000_fixed.gp
	gnuplot $<

$(OGCTGFIXEDGPS) : plots/$(OGID)_vs_%_ctgs_l100_fixed.gp : plots/$(OGID)_vs_%_ctgs_l100.gp
	../tools/mummerFixGP.awk -v "qid=$(OGID)" -v "rid=$*" $< > $@

$(OGFIXEDGPS) : plots/$(OGID)_vs_%_l100_fixed.gp : plots/$(OGID)_vs_%_l100.gp
	../tools/mummerFixGP.awk -v "qid=$(OGID)" -v "rid=$*" $< > $@

$(CTGFIXEDGPS) : plots/%_ctgs_vs_chroms_l100c1000_fixed.gp : plots/%_ctgs_vs_chroms_l100c1000.gp
	../tools/mummerFixGP.awk -v "qid=$* scaffolded contigs" -v "rid=$* arms" $< > $@

$(UNSCAFFIXEDGPS) : plots/%_unscaf_vs_chroms_l100c1000_fixed.gp : plots/%_unscaf_vs_chroms_l100c1000.gp
	../tools/mummerFixGP.awk -v "qid=$* unscaffolded contigs" -v "rid=$* arms" $< > $@

$(OLDFIXEDGPS) : plots/$(OLDID)_vs_%_l100c1000_fixed.gp : plots/$(OLDID)_vs_%_l100c1000.gp
	../tools/mummerFixGP.awk -v "qid=$(OLDID)" -v "rid=$*" $< > $@

$(MILLERFIXEDGPS) : plots/$(MILLER)_vs_%_l100c1000_fixed.gp : plots/$(MILLER)_vs_%_l100c1000.gp
	../tools/mummerFixGP.awk -v "qid=$(MILLER)" -v "rid=$*" $< > $@

$(CLOSEFIXEDGPS) : plots/$(CLOSEID)_vs_%_l100c1000_fixed.gp : plots/$(CLOSEID)_vs_%_l100c1000.gp
	../tools/mummerFixGP.awk -v "qid=$(CLOSEID)" -v "rid=$*" $< > $@

plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000_fixed.gp : plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.gp
	../tools/mummerFixGP.awk -v "qid=$(word 1,$(NEWASMS))" -v "rid=$(word 2,$(NEWASMS))" $< > $@

$(OGCTGGPS) : plots/$(OGID)_vs_%_ctgs_l100.gp : plot_order/$(OGID)_scaf_plot_order.tsv plot_order/%_contig_plot_order.tsv deltas/$(OGID)_vs_%_ctgs_l100.delta plots
	mummerplot --png --large --filter -Q $(word 1,$^) -R $(word 2,$^) -p plots/$(OGID)_vs_$*_ctgs_l100 $(word 3,$^) 2> logs/mummerplot_$(OGID)_vs_$*_ctgs_l100.stderr

$(OGGPS) : plots/$(OGID)_vs_%_l100.gp : plot_order/$(OGID)_scaf_plot_order.tsv plot_order/%_scaf_plot_order.tsv deltas/$(OGID)_vs_%_l100.delta plots
	mummerplot --png --large --filter -Q $(word 1,$^) -R $(word 2,$^) -p plots/$(OGID)_vs_$*_l100 $(word 3,$^) 2> logs/mummerplot_$(OGID)_vs_$*_l100.stderr

$(CTGGPS) : plots/%_ctgs_vs_chroms_l100c1000.gp : plot_order/%_contig_plot_order.tsv plot_order/%_scaf_plot_order.tsv deltas/%_ctgs_vs_chroms_l100c1000.delta plots
	mummerplot --png --large --filter -Q $(word 1,$^) -R $(word 2,$^) -p plots/$*_ctgs_vs_chroms_l100c1000 $(word 3,$^) 2> logs/mummerplot_$*_ctgs_vs_chroms_l100c1000.stderr

$(UNSCAFGPS) : plots/%_unscaf_vs_chroms_l100c1000.gp : plot_order/%_unscaf_plot_order.tsv plot_order/%_scaf_plot_order.tsv deltas/%_ctgs_vs_chroms_l100c1000.delta plots
	mummerplot --png --large --filter -Q $(word 1,$^) -R $(word 2,$^) -p plots/$*_unscaf_vs_chroms_l100c1000 $(word 3,$^) 2> logs/mummerplot_$*_unscaf_vs_chroms_l100c1000.stderr

$(OLDGPS) : plots/$(OLDID)_vs_%_l100c1000.gp : plot_order/$(OLDID)_scaf_plot_order.tsv plot_order/%_scaf_plot_order.tsv deltas/$(OLDID)_vs_%_l100c1000.delta plots
	mummerplot --png --large --filter -Q $(word 1,$^) -R $(word 2,$^) -p plots/$(OLDID)_vs_$*_l100c1000 $(word 3,$^) 2> logs/mummerplot_$(OLDID)_vs_$*_l100c1000.stderr

$(MILLERGPS) : plots/$(MILLER)_vs_%_l100c1000.gp : plot_order/$(MILLER)_scaf_plot_order.tsv plot_order/%_scaf_plot_order.tsv deltas/$(MILLER)_vs_%_l100c1000.delta plots
	mummerplot --png --large --filter -Q $(word 1,$^) -R $(word 2,$^) -p plots/$(MILLER)_vs_$*_l100c1000 $(word 3,$^) 2> logs/mummerplot_$(MILLER)_vs_$*_l100c1000.stderr

$(CLOSEGPS) : plots/$(CLOSEID)_vs_%_l100c1000.gp : plot_order/$(CLOSEID)_scaf_plot_order.tsv plot_order/%_scaf_plot_order.tsv deltas/$(CLOSEID)_vs_%_l100c1000.delta plots
	mummerplot --png --large --filter -Q $(word 1,$^) -R $(word 2,$^) -p plots/$(CLOSEID)_vs_$*_l100c1000 $(word 3,$^) 2> logs/mummerplot_$(CLOSEID)_vs_$*_l100c1000.stderr

plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.gp : plot_order/$(word 1,$(NEWASMS))_scaf_plot_order.tsv plot_order/$(word 2,$(NEWASMS))_scaf_plot_order.tsv deltas/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.delta plots
	mummerplot --png --large --filter -Q $(word 1,$^) -R $(word 2,$^) -p plots/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000 $(word 3,$^) 2> logs/mummerplot_$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.stderr

$(CTGORDERS) : plot_order/%_contig_plot_order.tsv : ../refs/contigs/%_fixed_configString.txt plot_order/%_contigs_filtered_softmasked_w60.fasta.fai
	../tools/contigPlotOrder.awk -v "chromlist=$(CHROMS)" $(word 1,$^) $(word 2,$^) > $@

$(UNSCAFORDERS) : plot_order/%_unscaf_plot_order.tsv : ../refs/contigs/%_fixed_configString.txt plot_order/%_contigs_filtered_softmasked_w60.fasta.fai
	../tools/contigPlotOrder.awk -v "nonchrom=1" -v "chromlist=$(CHROMS)" $(word 1,$^) $(word 2,$^) > $@

$(REFORDERS) : plot_order/%_scaf_plot_order.tsv : plot_order/%.fasta.fai
	../tools/mummerPlotOrder.awk -v "chromlist=$(CHROMS)" $< > $@

plot_order/%_contigs_filtered_softmasked_w60.fasta.fai : ../refs/contigs/%_contigs_filtered_softmasked_w60.fasta plot_order
	ln -sf ../$(word 1,$^) plot_order/$*_contigs_filtered_softmasked_w60.fasta; samtools faidx plot_order/$*_contigs_filtered_softmasked_w60.fasta

plot_order/$(OLDID).fasta.fai : ../refs/old/$(OLDID).fasta plot_order
	ln -sf ../$(word 1,$^) plot_order/$(OLDID).fasta; samtools faidx plot_order/$(OLDID).fasta

plot_order/$(MILLER).fasta.fai : ../refs/old/$(MILLER).fasta plot_order
	ln -sf ../$(word 1,$^) plot_order/$(MILLER).fasta; samtools faidx plot_order/$(MILLER).fasta

plot_order/%.fasta.fai : ../refs/%.fasta plot_order
	ln -sf ../$(word 1,$^) plot_order/$*.fasta; samtools faidx plot_order/$*.fasta

$(OGCTGDELTAS) : deltas/$(OGID)_vs_%_ctgs_l100.delta : ../refs/contigs/%_contigs_filtered_softmasked_w60.fasta ../refs/$(OGID).fasta logs deltas
	/usr/bin/time -v nucmer -l100 -p deltas/$(OGID)_vs_$*_ctgs_l100 $(word 1,$^) $(word 2,$^) 2> logs/nucmer_$(OGID)_vs_$*_ctgs_l100.stderr > logs/nucmer_$(OGID)_vs_$*_ctgs_l100.stdout

$(OGDELTAS) : deltas/$(OGID)_vs_%_l100.delta : ../refs/%.fasta ../refs/$(OGID).fasta logs deltas
	/usr/bin/time -v nucmer -l100 -p deltas/$(OGID)_vs_$*_l100 $(word 1,$^) $(word 2,$^) 2> logs/nucmer_$(OGID)_vs_$*_l100.stderr > logs/nucmer_$(OGID)_vs_$*_l100.stdout

$(CTGDELTAS) : deltas/%_ctgs_vs_chroms_l100c1000.delta : ../refs/%.fasta ../refs/contigs/%_contigs_filtered_softmasked_w60.fasta logs deltas
	/usr/bin/time -v nucmer -l100 -c1000 -p deltas/$*_ctgs_vs_chroms_l100c1000 $(word 1,$^) $(word 2,$^) 2> logs/nucmer_$*_ctgs_vs_chroms_l100c1000.stderr > logs/nucmer_$*_ctgs_vs_chroms_l100c1000.stdout

$(OLDDELTAS) : deltas/$(OLDID)_vs_%_l100c1000.delta : ../refs/%.fasta ../refs/old/$(OLDID).fasta logs deltas
	/usr/bin/time -v nucmer -l100 -c1000 -p deltas/$(OLDID)_vs_$*_l100c1000 $(word 1,$^) $(word 2,$^) 2> logs/nucmer_$(OLDID)_vs_$*_l100c1000.stderr > logs/nucmer_$(OLDID)_vs_$*_l100c1000.stdout

deltas/$(MILLER)_vs_$(OLDID)_l100c1000.delta : ../refs/old/$(OLDID).fasta ../refs/old/$(MILLER).fasta logs deltas
	/usr/bin/time -v nucmer -l100 -c1000 -p deltas/$(MILLER)_vs_$(OLDID)_l100c1000 $(word 1,$^) $(word 2,$^) 2> logs/nucmer_$(MILLER)_vs_$(OLDID)_l100c1000.stderr > logs/nucmer_$(MILLER)_vs_$(OLDID)_l100c1000.stdout

deltas/$(MILLER)_vs_$(MILLERMATCH)_l100c1000.delta : ../refs/$(MILLERMATCH).fasta ../refs/old/$(MILLER).fasta logs deltas
	/usr/bin/time -v nucmer -l100 -c1000 -p deltas/$(MILLER)_vs_$(MILLERMATCH)_l100c1000 $(word 1,$^) $(word 2,$^) 2> logs/nucmer_$(MILLER)_vs_$(MILLERMATCH)_l100c1000.stderr > logs/nucmer_$(MILLER)_vs_$(MILLERMATCH)_l100c1000.stdout

$(CLOSEDELTAS) : deltas/$(CLOSEID)_vs_%_l100c1000.delta : ../refs/%.fasta ../refs/$(CLOSEID).fasta logs deltas
	/usr/bin/time -v nucmer -l100 -c1000 -p deltas/$(CLOSEID)_vs_$*_l100c1000 $(word 1,$^) $(word 2,$^) 2> logs/nucmer_$(CLOSEID)_vs_$*_l100c1000.stderr > logs/nucmer_$(CLOSEID)_vs_$*_l100c1000.stdout

deltas/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.delta : ../refs/$(word 1,$(NEWASMS)).fasta ../refs/$(word 2,$(NEWASMS)).fasta logs deltas
	/usr/bin/time -v nucmer -l100 -c1000 -p deltas/$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000 $(word 2,$^) $(word 1,$^) 2> logs/nucmer_$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.stderr > logs/nucmer_$(word 1,$(NEWASMS))_vs_$(word 2,$(NEWASMS))_l100c1000.stdout

$(SUBDIRS) :
	mkdir -p $@

clean :
	for i in $(SUBDIRS); do rm -f $${i}/*.fasta $${i}/*.fai $${i}/*.tsv $${i}/*.delta $${i}/*.[fr]plot $${i}/*.gp $${i}/*.filter $${i}/*.png $${i}/*.stderr $${i}/*.stdout; [[ ! -d $${i} ]] || rmdir $${i}; done
