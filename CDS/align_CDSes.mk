SHELL=/bin/bash

#Align CDSes (reference CDSes with PRANK, then add pseudoref CDSes)
#Assumes refs and annotations are available in subdirectories
# of parent directory, plus pseudoreferences

#Expects the following in your PATH:
#prank
#fasta_formatter
#The rest of these should be accessible via relative paths in the tools
# subdirectory:
#addAlignedGaps.pl
#constructCDSesFromGFF3.pl
#fakeHaplotype.pl
#parseFASTARecords.pl
#omitInbredHaplotype.awk
#trimAlignment.awk
#trimCDSes.awk
#trimSCOmapHeaders.awk
#These Perl and awk scripts are available in the tools subdirectory
# of the Github repository.

#Options to adjust:
#Map of single copy orthogroups:
#(prep_proteomes.mk creates a file called Orthogroups_SingleCopy_renamed.tsv
# in its directory (i.e. where you run OrthoFinder, though not OrthoFinder's
# working directory))
SCOMAP := ../OrthoFinder_clustalo/Orthogroups_SingleCopy_renamed.tsv
#References to exclude from alignments:
EXCLUDEDREFS := Dyak_Tai18E2

#DO NOT CHANGE THE BELOW:
#Temporary single copy orthogroup map for use with pseudorefs:
SAMPLESCOMAP := Orthogroups_SingleCopy_renamed_forPseudorefs.tsv

#Subdirectories that need to be created:
SUBDIRS := unparsed_CDSes_refs unparsed_CDSes_samples raw_CDSes_refs raw_CDSes_samples trimmed_unaligned_refs trimmed_unaligned_samples trimmed_aligned_refs trimmed_aligned_all

#References inferred from directory above the current:
REFS := $(filter-out $(EXCLUDEDREFS),$(basename $(notdir $(wildcard ../refs/*.fasta))))
CDSREFS := $(addprefix unparsed_CDSes_refs/,$(addsuffix .fasta,$(REFS)))

#Run this part after parsed_CDSes_refs is complete:
PSEUDOREFS := $(basename $(notdir $(wildcard ../pseudorefs/*.fasta)))
CDSSAMPLES := $(addprefix unparsed_CDSes_samples/,$(addsuffix .fasta,$(PSEUDOREFS)))
OGs := $(patsubst %_unwrapped.fasta,%,$(notdir $(wildcard raw_CDSes_refs/OG*.fasta)))
RAWREFS := $(addprefix raw_CDSes_refs/,$(addsuffix _unwrapped.fasta,$(OGs)))
RAWSAMPLES := $(addprefix raw_CDSes_samples/,$(addsuffix _unwrapped.fasta,$(OGs)))
TRIMMEDREFS := $(addprefix trimmed_unaligned_refs/,$(addsuffix _trimmed.fasta,$(OGs)))
TRIMMEDSAMPLES := $(addprefix trimmed_unaligned_samples/,$(addsuffix _trimmed.fasta,$(OGs)))
ALIGNEDREFS := $(addprefix trimmed_aligned_refs/,$(addsuffix _prank.best.fas,$(OGs)))
ALIGNEDALL := $(addprefix trimmed_aligned_all/,$(addsuffix _trimmed.fasta,$(OGs)))

.PHONY : check clean usage

.SECONDARY : $(SAMPLESCOMAP) $(CDSREFS) $(CDSSAMPLES) $(RAWREFS) $(RAWSAMPLES) $(TRIMMEDREFS) $(TRIMMEDSAMPLES) $(ALIGNEDREFS) $(ALIGNEDALL)

usage :
	@echo "Usage:"
	@echo "make -f align_CDSes.mk [task]"
	@echo "Tasks:"
	@echo "parse_ref_CDSes -> Extract CDSes from refs, and parse into Orthogroups"
	@echo "align_refs -> (Run after parse_ref_CDSes) Trim ref CDSes and align with PRANK"
	@echo "parse_sample_CDSes -> (Run after align_refs) Extract CDSes from pseudorefs, and parse into Orthogroups"
	@echo "align_all -> (Run after parse_sample_CDSes) Merge pseudoref haplotypes into alignment and trim"
	@echo "check -> Output diagnostic TSVs of the number of unaligned haplotypes per Orthogroup, and the number of haplotypes per Orthogroup after merging alignments and trimming"
	@echo "clean -> Clean up all output files and directories"

check : samples_per_OG_unaligned.tsv trimmed_aligned_all_number_of_seqs.tsv

#Check number of pseudorefs per OG:
samples_per_OG_unaligned.tsv : parse_sample_CDSes
	for i in $(OGs); do printf "$${i}\t"; fgrep -c ">" raw_CDSes_samples/$${i}.fasta; done > $@

#Output number of rows in final alignments:
trimmed_aligned_all_number_of_seqs.tsv : align_all
	for i in $(ALIGNEDALL); do fn=`basename $${i}`; orthogroup=$${fn%_trimmed.fasta}; printf "$${orthogroup}\t"; fgrep -c ">" $${i}; done > $@

#Node in DAG parent to all aligned CDSes:
#Note: We don't use make auto variable here as it would give a
# bash arg list error (execvp: /bin/bash: Argument list too long)
align_all : $(ALIGNEDALL)
	@echo "Done parsing, trimming, and aligning pseudoref CDSes, and trimming the overall alignments"
	fgrep -c ">" trimmed_aligned_all/*_trimmed.fasta | cut -d":" -f2 | sort | uniq -c > $@

#Add to ref aligned CDSes, and trim alignments:
trimmed_aligned_all/%_trimmed.fasta : trimmed_aligned_refs/%_prank.best.fas trimmed_unaligned_samples/%_trimmed.fasta trimmed_aligned_all
	../tools/addAlignedGaps.pl -a <(fasta_formatter -i $(word 1,$^)) -i $(word 2,$^) | ../tools/trimAlignment.awk > $@

#Split CDSes into fake haplotypes and trim:
trimmed_unaligned_samples/OG%_trimmed.fasta : raw_CDSes_samples/OG%_unwrapped.fasta trimmed_unaligned_samples
	cat <(fasta_formatter -i $(word 1,$^) | ../tools/fakeHaplotype.pl -b) <(fasta_formatter -i $(word 1,$^) | ../tools/fakeHaplotype.pl -b -a | ../tools/omitInbredHaplotype.awk) | ../tools/trimCDSes.awk > $@

#Parse CDSes by OG:
parse_sample_CDSes : $(SAMPLESCOMAP) $(CDSSAMPLES) raw_CDSes_samples
	@echo "Parsing pseudoref CDSes by orthogroup"
	cd raw_CDSes_samples; \
	mkdir -p logs; \
	ls $(addprefix ../,$(CDSSAMPLES)) > unparsed_CDSes_samples.fofn; \
	/usr/bin/time -v ../../tools/parseFASTARecords.pl -f unparsed_CDSes_samples.fofn -m $(addprefix ../,$(word 1,$^)) -dddd 2> logs/parseFASTA.stderr > logs/parseFASTA.stdout; \
	cd ..; \
	fgrep -c ">" raw_CDSes_samples/OG* | cut -d":" -f2 | sort | uniq -c > $@

#Extract CDSes from pseudorefs:
unparsed_CDSes_samples/%.fasta : ../pseudorefs/%.fasta unparsed_CDSes_samples
	mkdir -p logs; \
	species=$$(echo "$*" | cut -d"_" -f1); \
	../tools/constructCDSesFromGFF3.pl -l -i $(word 1,$^) -g ../annotations/$${species}.gff3 -p $* > $@ 2> logs/txome_$*.stderr

#Node in DAG parent to all aligned reference CDSes:
#Note: We don't use make auto variable here as it would give a
# bash arg list error (execvp: /bin/bash: Argument list too long)
align_refs : $(ALIGNEDREFS)
	@echo "Done parsing, trimming, and aligning SCO CDSes of refs"
	fgrep -c ">" trimmed_aligned_refs/*_prank.best.fas | cut -d":" -f2 | sort | uniq -c > $@

#Align ref CDSes with PRANK:
trimmed_aligned_refs/%_prank.best.fas : trimmed_unaligned_refs/%_trimmed.fasta trimmed_aligned_refs
	mkdir -p prank_logs; \
	/usr/bin/time -v prank -d=$(word 1,$^) -o=trimmed_aligned_refs/$*_prank -codon -F -verbose 2> prank_logs/$*_prank.stderr > prank_logs/$*_prank.stdout

#Trim ref CDSes:
trimmed_unaligned_refs/OG%_trimmed.fasta : raw_CDSes_refs/OG%_unwrapped.fasta trimmed_unaligned_refs
	../tools/trimCDSes.awk $(word 1,$^) > $@

#Parse CDSes by OG:
parse_ref_CDSes : $(SAMPLESCOMAP) $(CDSREFS) raw_CDSes_refs
	@echo "Parsing ref CDSes by orthogroup"
	cd raw_CDSes_refs; \
	mkdir -p logs; \
	ls $(addprefix ../,$(CDSREFS)) > unparsed_CDSes_refs.fofn; \
	/usr/bin/time -v ../../tools/parseFASTARecords.pl -f unparsed_CDSes_refs.fofn -m $(addprefix ../,$(word 1,$^)) -dddd 2> logs/parseFASTA.stderr > logs/parseFASTA.stdout; \
	cd ..; \
	fgrep -c ">" raw_CDSes_refs/OG* | cut -d":" -f2 | sort | uniq -c > $@

#Trim SCO map headers and excluded requested refs:
$(SAMPLESCOMAP) : $(SCOMAP)
	@echo "Trimming single copy orthogroup map headers for use with pseudorefs"
	../tools/trimSCOmapHeaders.awk -v "excludedcols=$(EXCLUDEDREFS)" $< > $@

#Extract CDSes from refs:
unparsed_CDSes_refs/%.fasta : ../refs/%.fasta ../annotations/%.gff3 unparsed_CDSes_refs
	mkdir -p logs; \
	/usr/bin/time -v ../tools/constructCDSesFromGFF3.pl -l -i $(word 1,$^) -g $(word 2,$^) -p $* > $@ 2> logs/txome_$*.stderr

$(SUBDIRS) :
	mkdir -p $@

clean :
	for i in $(SUBDIRS); do rm -f $${i}/*.fasta $${i}/logs/* $${i}/*.fofn $${i}/*.fas; [[ ! -d $${i}/logs ]] || rmdir $${i}/logs/; [[ ! -d $${i} ]] || rmdir $${i}; done
	rm -f prank_logs/* logs/*
	[[ ! -d prank_logs ]] || rmdir prank_logs
	[[ ! -d logs ]] || rmdir logs
	rm -f $(SAMPLESCOMAP)
	rm -f parse_ref_CDSes align_refs parse_sample_CDSes align_all
	rm -f samples_per_OG_unaligned.tsv trimmed_aligned_all_number_of_seqs.tsv
