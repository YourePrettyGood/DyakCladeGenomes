SHELL=/bin/bash

#Generate length-ordered contig length distributions for
# all assemblies, and compare to Dmel scaffold lengths.
#Plot generated at the end is a cumulative contig length
# distribution.
#Assumes refs are available in subdirectory of parent directory

#Expects the following in your PATH:
#fasta_formatter
#The rest of these should be accessible via relative paths in the tools
# subdirectory:
#AssemblyContiguityPlot.R
#This R script is available in the tools subdirectory
# of the Github repository.

#Options to adjust:
#Ref to use scaffolds from as ideal:
IDEALREF := Dmel_ISO1
#Final plot filename:
PLOTNAME := MSTYY_ctg_len_dist.pdf
#Exclude a ref:
EXCLUDEREF := Dsim_w501

#DO NOT CHANGE THE BELOW:
#References inferred from directory above the current:
REFS := $(filter-out $(EXCLUDEREF),$(basename $(notdir $(wildcard ../refs/*.fasta))))
#Old reference(s) (not PacBio one(s)), inferred from "old" subdirectory of
# "refs" directory above current:
OLDREFS := $(basename $(notdir $(wildcard ../refs/old/*.fasta)))
#Contig length files:
CTGLENS := $(addsuffix _ctglens.tsv,$(REFS))
#Scaffold length file:
SCAFLEN := $(addsuffix _scaflens.tsv,$(IDEALREF))
#Contig length file(s) for old references (not our PacBio ones):
EXTRACTGLENS := $(addsuffix _ctglens.tsv,$(OLDREFS))

.PHONY : all clean usage

.SECONDARY : $(CTGLENS) $(SCAFLEN) $(EXTRACTGLENS) combined_lens.tsv

usage :
	@echo "Usage:"
	@echo "make -f ctg_lengths.mk [task]"
	@echo "Tasks:"
	@echo "all -> Provide length-sorted lists of contig lengths (plus scaffold lengths for an ideal genome)"
	@echo "clean -> Clean up all output files and directories"

all : $(PLOTNAME)

#Generate the cumulative contig length plot:
$(PLOTNAME) : combined_lens.tsv
	../tools/AssemblyContiguityPlot.R $@ $< Ideal_$(IDEALREF) $(REFS) $(OLDREFS)

#Combine the contig and scaffold length files:
combined_lens.tsv : $(CTGLENS) $(SCAFLEN) $(EXTRACTGLENS)
	cat $^ > $@

$(EXTRACTGLENS) : %_ctglens.tsv : ../refs/old/%.fasta
	fasta_formatter -i $< | awk '/^>/{header=$$0;}!/^>/{split($$0, ctgs, /[N]+/); for (ctg in ctgs) {print header"_ctg"ctg; print ctgs[ctg];};}' | awk '!/^>/{print length($$0);}' | sort -k1,1nr | awk '{print "$*\t"NR"\t"$$1;}' > $@

#Sort scaffold lengths for the ideal reference:
$(SCAFLEN) : ../refs/$(IDEALREF).fasta
	fasta_formatter -i $< | awk '!/^>/{print length($$0);}' | sort -k1,1nr | awk '{print "Ideal_$(IDEALREF)\t"NR"\t"$$1;}' > $@

#Sort contig lengths for the other references:
$(CTGLENS) : %_ctglens.tsv : ../refs/%.fasta
	fasta_formatter -i $< | awk '/^>/{header=$$0;}!/^>/{split($$0, ctgs, /[N]+/); for (ctg in ctgs) {print header"_ctg"ctg; print ctgs[ctg];};}' | awk '!/^>/{print length($$0);}' | sort -k1,1nr | awk '{print "$*\t"NR"\t"$$1;}' > $@

clean :
	rm -f $(CTGLENS) $(SCAFLEN) $(EXTRACTGLENS) combined_lens.tsv $(PLOTNAME)
