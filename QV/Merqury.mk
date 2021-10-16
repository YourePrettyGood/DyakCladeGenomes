SHELL=/bin/bash

#Assess per-base accuracy of the assemblies using Merqury.
#Assumes refs and asm_reads are available in subdirectory
# of parent directory

#Expects the following to be set up in your environment:
# meryl v1.1 in your PATH
# MERQURY environment variable set to the path to the Merqury install
#

#Options to adjust:
#Final QV summary filename:
QVSUMMARY := DyakClade_Merqury_QVs.tsv
#Final completeness summary filename:
COMPLETENESS := DyakClade_Merqury_completeness.tsv
#Exclude refs:
EXCLUDEREF := Dmel_ISO1 Dsim_w501
#Size of k-mers (selected with best_k.sh, but we're hard-coding for now:
K := 19
#Memory for meryl to use (in GB):
MEM := 32
#Number of threads for meryl to use:
NTHREADS := 16
#Reads to use for "old" assemblies:
OLDREADSET := Dyak_Tai18E2

#DO NOT CHANGE THE BELOW:
#References inferred from directory above the current:
REFS := $(filter-out $(EXCLUDEREF),$(basename $(notdir $(wildcard ../refs/*.fasta))))
#Old reference(s) (not PacBio one(s)), inferred from "old" subdirectory of
# "refs" directory above current:
OLDREFS := $(basename $(notdir $(wildcard ../refs/old/*.fasta)))
#
QVS := $(addsuffix _qv.tsv,$(REFS))
OLDQVS := $(addsuffix _qv.tsv,$(OLDREFS))
COMPLETENESSES := $(addsuffix _completeness.tsv,$(REFS))
OLDCOMPLETENESSES := $(addsuffix _completeness.tsv,$(OLDREFS))
READDBS := $(addprefix DBs/reads/,$(addsuffix .meryl,$(REFS)))
OLDREADDBS := $(addprefix DBs/reads/,$(addsuffix .meryl,$(OLDREFS)))
#Subdirectories to create:
SUBIDRS := DBs/asm DBs/perFASTQ DBs/reads

.PHONY : all clean usage

.SECONDARY :

usage :
	@echo "Usage:"
	@echo "make -f ctg_lengths.mk [task]"
	@echo "Tasks:"
	@echo "all -> Calculate QV and kmer completeness for Dyak clade assemblies"
	@echo "clean -> Clean up all output files and directories"

all : $(QVSUMMARY) $(COMPLETENESS)

$(QVSUMMARY) : $(QVS) $(OLDQVS)
	printf "Assembly\tAssembly-only kmers\tTotal kmers in assembly\tQV\tProbability of error\n" > $@; \
	cat $^ >> $@

$(QVS) $(OLDQVS) : %_qv.tsv : DBs/asms/%.0.meryl DBs/asms/%.meryl
	../tools/merquryQV.awk -v "k=$(K)" -v "asm=$*" <(meryl statistics $(word 1,$^)) <(meryl statistics $(word 2,$^)) > $@

$(COMPLETENESS) : $(COMPLETENESSES) $(OLDCOMPLETENESSES)
	printf "Assembly\tReliable read kmers in assembly\tTotal reliable kmers\tKmer completeness\n" > $@; \
	cat $^ >> $@

$(COMPLETENESSES) $(OLDCOMPLETENESSES) : %_completeness.tsv : DBs/asm/%.solid.meryl DBs/reads/%.meryl
	../tools/merquryCompleteness.awk -v "asm=$*" <(meryl statistics $(word 1,$^)) <(meryl statistics $(word 2,$^)) > $@

DBs/asm/%.0.meryl : DBs/asm/%.meryl DBs/reads/%.meryl
	meryl difference output $@ $(word 1,$^) $(word 2,$^)

DBs/asm/%.solid.meryl : DBs/asm/%.meryl DBs/reads/%.reliable.meryl
	meryl intersect output $@ $(word 1,$^) $(word 2,$^)

DBs/reads/%.reliable.meryl : DBs/reads/%.meryl DBs/reads/%.hist.ploidy
	meryl greater-than $$(awk 'NR==2{print $$NF;}' $(word 2,$^)) output $@ $(word 1,$^)

DBs/reads/%.hist.ploidy : DBs/reads/%.hist
	java -Xmx1g -jar ${MERQURY}/eval/kmerHistToPloidyDepth.jar $< > $@

DBs/reads/%.hist : DBs/reads/%.meryl
	meryl histogram $< > $@

DBs/asm/%.meryl : ../refs/%.fasta DBs/asm
	meryl k=$(K) threads=$(NTHREADS) memory=$(MEM) count output $@ $<

$(OLDREADDBS) : DBs/reads/%.meryl : DBs/reads/$(OLDREADSET).meryl
	ln -sf $(notdir $@) $<

$(READDBS) : DBs/reads/%.meryl : ../asm_reads/% DBs/reads DBs/perFASTQ
	for f in $</*.fastq.gz; do \
	   fn=$$(basename $${f}); \
	   prefix=$${fn%.fastq.gz}; \
	   meryl k=$(K) memory=$(MEM) threads=$(NTHREADS) count $${f} output $(word 3,$^)/$${prefix}.meryl; \
	done; \
	meryl k=$(K) memory=$(MEM) threads=$(NTHREADS) union-sum output $(word 2,$^).meryl $(word 1,$^)/*.meryl

$(SUBDIRS) : 
	mkdir -p $@

#for asm in Dsan_STOCAGO1482 Dtei_GT53w Dyak_NY73PB Dyak_Tai18E2 Dyak_caf1;
#   do
#   if [[ "${asm}" == "Dyak_caf1" ]]; then
#      asmfn="../../refs/old/${asm}.fasta";
#      reads="../DBs/reads/reads_Dyak_Tai18E2.meryl";
#   else
#      asmfn="../../refs/${asm}.fasta";
#      reads="../DBs/reads/reads_${asm}.meryl";
#   fi;
#
#done

clean :
	rm -f $(COMPLETENESS) $(COMPLETENESSES) $(OLDCOMPLETENESS) $(QVSUMMARY) $(QVS) $(OLDQVS)
	for i in $(SUBDIRS); do for d in $${i}/*; do rm -f $${d}/*.merylData $${d}/*.merylIndex $${d}/merylIndex; [[ ! -d "$${d}" ]] && rmdir $${d}; done; [[ ! -d "$${i}" ]] && rmdir $${i}; done
