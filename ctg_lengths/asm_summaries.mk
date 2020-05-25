SHELL=/bin/bash

#Generate assembly summary statistics for Table 1
#Assumes refs are available in subdirectory of parent directory

#Expects the following in your PATH:
#fasta_formatter
#The rest of these should be accessible via relative paths in the tools
# subdirectory:
#FASTAtoContigs.awk
#NX.pl
#These awk and perl scripts are available in the tools subdirectory
# of the Github repository.

#Options to adjust:
SUMMARYFILE := DyakClade_asm_summary_stats.tsv

#DO NOT CHANGE THE BELOW:
#References inferred from directory above the current:
REFS := $(basename $(notdir $(wildcard ../refs/*.fasta)))
#Old reference(s) (not PacBio one(s)), inferred from "old" subdirectory of
# "refs" directory above current:
OLDREFS := $(basename $(notdir $(wildcard ../refs/old/*.fasta)))
#Assembly summary statistics files:
SUMMARIES := $(addsuffix _asm_stats.tsv,$(REFS) $(OLDREFS))
#Scaffold summary statistics:
REFSCAFSUMMARIES := $(addsuffix _scaf_stats.tsv,$(REFS))
OLDREFSCAFSUMMARIES := $(addsuffix _scaf_stats.tsv,$(OLDREFS))
#Contig summary statistics:
REFCTGSUMMARIES := $(addsuffix _ctg_stats.tsv,$(REFS))
OLDREFCTGSUMMARIES := $(addsuffix _ctg_stats.tsv,$(OLDREFS))

.PHONY : all distclean clean usage

usage :
	@echo "Usage:"
	@echo "make -f asm_summaries.mk [task]"
	@echo "Tasks:"
	@echo "all -> Calculate assembly summary statistics"
	@echo "distclean -> Clean up intermediate files"
	@echo "clean -> Clean up all output files and directories"

all : $(SUMMARYFILE)

#Combine the individual assembly summary files:
$(SUMMARYFILE) : $(SUMMARIES)
	cat <(printf "Asm\tAsmSize\tCtgN50\tCtgN90\tCtgL50\tCtgL90\tNumCtgs\tScafN50\tScafN90\tScafL50\tScafL90\tNumScafs\n") $^ > $@

#Combine scaffold and contig summaries:
$(SUMMARIES) : %_asm_stats.tsv : %_ctg_stats.tsv %_scaf_stats.tsv
	paste <(printf "$*\n") $(word 1,$^) $(word 2,$^) > $@

#Calculate scaffold summaries for refs:
$(REFSCAFSUMMARIES) : %_scaf_stats.tsv : ../refs/%.fasta
	fasta_formatter -i $< | ../tools/NX.pl -t - 50,90 | cut -f2,3,4,5,8 | tail -n1 > $@

#Calculate scaffold summaries for old refs:
$(OLDREFSCAFSUMMARIES) : %_scaf_stats.tsv : ../refs/old/%.fasta
	fasta_formatter -i $< | ../tools/NX.pl -t - 50,90 | cut -f2,3,4,5,8 | tail -n1 > $@

#Calculate contig summaries for refs:
$(REFCTGSUMMARIES) : %_ctg_stats.tsv : ../refs/%.fasta
	fasta_formatter -i $< | ../tools/FASTAtoContigs.awk | ../tools/NX.pl -t - 50,90 | cut -f1,2,3,4,5,8 | tail -n1 > $@

#Calculate contig summaries for old refs:
$(OLDREFCTGSUMMARIES) : %_ctg_stats.tsv : ../refs/old/%.fasta
	fasta_formatter -i $< | ../tools/FASTAtoContigs.awk | ../tools/NX.pl -t - 50,90 | cut -f1,2,3,4,5,8 | tail -n1 > $@

distclean :
	rm -f $(SUMMARIES) $(REFSCAFSUMMARIES) $(OLDREFSCAFSUMMARIES) $(REFCTGSUMMARIES) $(OLDREFCTGSUMMARIES)

clean :
	rm -f $(SUMMARYFILE) $(SUMMARIES) $(REFCTGSUMMARIES) $(OLDREFCTGSUMMARIES) $(REFSCAFSUMMARIES) $(OLDREFSCAFSUMMARIES)
