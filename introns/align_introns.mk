SHELL=/bin/bash

#Align introns (reference introns with Clustal Omega, then add
# pseudoref introns)
#Assumes refs and annotations are available in subdirectories
# of parent directory, plus pseudoreferences

#Expects the following in your PATH:
#clustalo
#fasta_formatter
#The rest of these should be accessible via relative paths in the tools
# subdirectory:
#addAlignedGaps.pl
#extractIntronsFromGFF3.pl
#fakeHaplotype.pl
#findSingleCopyOrthologIntrons.pl
#parseFASTARecords.pl
#omitInbredHaplotype.awk
#trimAlignment.awk
#trimIntrons.awk
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
#Maximum relative difference in flanking exon lengths to accept
# for assigning orthology of introns:
THRESHOLD := 0.20

#DO NOT CHANGE THE BELOW:
#Temporary single copy orthogroup map for use with pseudorefs:
SAMPLESCOMAP := Orthogroups_SingleCopy_renamed_forPseudorefs.tsv
#Intron orthology map:
INTRONMAP := Orthologous_introns.tsv
#Temporary intron orthology map for use with pseudorefs:
SAMPLEINTRONMAP := Orthologous_introns_forPseudorefs.tsv

#Subdirectories that need to be created:
SUBDIRS := unparsed_introns_refs unparsed_introns_samples raw_introns_refs raw_introns_samples trimmed_unaligned_refs trimmed_unaligned_samples trimmed_aligned_refs trimmed_aligned_all

#References inferred from directory above the current:
REFS := $(filter-out $(EXCLUDEDREFS),$(basename $(notdir $(wildcard ../refs/*.fasta))))
INTRONREFS := $(addprefix unparsed_introns_refs/,$(addsuffix .fasta,$(REFS)))

#Run this part after parsed_introns_refs is complete:
PSEUDOREFS := $(basename $(notdir $(wildcard ../pseudorefs/*.fasta)))
INTRONSAMPLES := $(addprefix unparsed_introns_samples/,$(addsuffix .fasta,$(PSEUDOREFS)))
introns := $(patsubst %_unwrapped.fasta,%,$(notdir $(wildcard raw_introns_refs/OG*.fasta)))
RAWREFS := $(addprefix raw_introns_refs/,$(addsuffix _unwrapped.fasta,$(introns)))
RAWSAMPLES := $(addprefix raw_introns_samples/,$(addsuffix _unwrapped.fasta,$(introns)))
TRIMMEDREFS := $(addprefix trimmed_unaligned_refs/,$(addsuffix _trimmed.fasta,$(introns)))
TRIMMEDSAMPLES := $(addprefix trimmed_unaligned_samples/,$(addsuffix _trimmed.fasta,$(introns)))
ALIGNEDREFS := $(addprefix trimmed_aligned_refs/,$(addsuffix _aligned.fasta,$(introns)))
ALIGNEDALL := $(addprefix trimmed_aligned_all/,$(addsuffix _trimmed.fasta,$(introns)))

.PHONY : check clean usage

.SECONDARY : $(SAMPLESCOMAP) $(INTRONMAP) $(SAMPLEINTRONMAP) $(INTRONREFS) $(INTRONSAMPLES) $(RAWREFS) $(RAWSAMPLES) $(TRIMMEDREFS) $(TRIMMEDSAMPLES) $(ALIGNEDREFS) $(ALIGNEDALL)

usage :
	@echo "Usage:"
	@echo "make -f align_introns.mk [task]"
	@echo "Tasks:"
	@echo "parse_ref_introns -> Extract introns from refs, and parse into orthologous intron groups"
	@echo "align_refs -> (Run after parse_ref_introns) Trim ref introns and align with Clustal Omega"
	@echo "parse_sample_introns -> (Run after align_refs) Extract introns from pseudorefs, and parse into orthologous intron groups"
	@echo "align_all -> (Run after parse_sample_introns) Merge pseudoref haplotypes into alignment and trim"
	@echo "check -> Output diagnostic TSVs of the number of unaligned haplotypes per intron group, and the number of haplotypes per intron group after merging alignments and trimming"
	@echo "clean -> Clean up all output files and directories"

check : samples_per_OIG_unaligned.tsv trimmed_aligned_all_number_of_seqs.tsv

#Check number of pseudorefs per orthologous intron group:
samples_per_OIG_unaligned.tsv : parsed_introns_samples
	for i in $(introns); do printf "$${i}\t"; fgrep -c ">" raw_introns_samples/$${i}.fasta; done > $@

#Output number of rows in final alignments:
trimmed_aligned_all_number_of_seqs.tsv : align_all
	for i in $(ALIGNEDALL); do fn=`basename $${i}`; introngroup=$${fn%_trimmed.fasta}; printf "$${introngroup}\t"; fgrep -c ">" $${i}; done > $@

align_all : $(ALIGNEDALL)
	@echo "Done parsing, trimming, and aligning pseudoref introns, and trimming the overall alignments"
	fgrep -c ">" trimmed_aligned_all/*_trimmed.fasta | cut -d":" -f2 | sort | uniq -c > $@

#Add to ref aligned introns, and trim alignments:
trimmed_aligned_all/%_trimmed.fasta : trimmed_aligned_refs/%_aligned.fasta trimmed_unaligned_samples/%_trimmed.fasta trimmed_aligned_all
	../tools/addAlignedGaps.pl -a <(fasta_formatter -i $(word 1,$^)) -i $(word 2,$^) 2> logs/addAlignedGaps_$*.stderr | ../tools/trimAlignment.awk > $@

#Split introns into fake haplotypes and trim:
trimmed_unaligned_samples/OG%_trimmed.fasta : raw_introns_samples/OG%_unwrapped.fasta trimmed_unaligned_samples
	cat <(fasta_formatter -i $(word 1,$^) | ../tools/fakeHaplotype.pl -b) <(fasta_formatter -i $(word 1,$^) | ../tools/fakeHaplotype.pl -b -a | ../tools/omitInbredHaplotype.awk) | ../tools/trimIntrons.awk > $@

#Parse introns by orthologous intron group:
parse_sample_introns : $(SAMPLEINTRONMAP) $(INTRONSAMPLES) raw_introns_samples
	@echo "Parsing pseudoref introns by orthologous intron group"
	cd raw_introns_samples; \
	mkdir -p logs; \
	ls $(addprefix ../,$(INTRONSAMPLES)) > unparsed_introns_samples.fofn; \
	/usr/bin/time -v ../../tools/parseFASTARecords.pl -f unparsed_introns_samples.fofn -m ../$(word 1,$^) -dddd 2> logs/parseFASTA.stderr > logs/parseFASTA.stdout; \
	cd ..; \
	fgrep -c ">" raw_introns_samples/*_unwrapped.fasta | cut -d":" -f2 | sort | uniq -c > $@

#Extract introns from pseudorefs:
unparsed_introns_samples/%.fasta : ../pseudorefs/%.fasta unparsed_introns_samples
	mkdir -p logs; \
	species=$$(echo "$*" | cut -d"_" -f1); \
	../tools/extractIntronsFromGFF3.pl -d -b CDS -i $(word 1,$^) -g ../annotations/$${species}.gff3 -p $* 2> logs/intronome_$*.stderr > $@

#Node in DAG parent to all aligned reference introns:
align_refs : $(ALIGNEDREFS)
	@echo "Done parsing, trimming, and aligning orthologous introns of refs"
	fgrep -c ">" trimmed_aligned_refs/*_aligned.fasta | cut -d":" -f2 | sort | uniq -c > $@

#Align ref introns with Clustal Omega:
trimmed_aligned_refs/%_aligned.fasta : trimmed_unaligned_refs/%_trimmed.fasta trimmed_aligned_refs
	mkdir -p clustalo_logs; \
	/usr/bin/time -v clustalo -i $(word 1,$^) 2> clustalo_logs/$*_clustalo.stderr | fasta_formatter -o $@

#Trim ref introns:
trimmed_unaligned_refs/OG%_trimmed.fasta : raw_introns_refs/OG%_unwrapped.fasta trimmed_unaligned_refs
	../tools/trimIntrons.awk $(word 1,$^) > $@

#Parse introns by orthologous intron group:
parse_ref_introns : $(INTRONMAP) $(INTRONREFS) raw_introns_refs
	@echo "Parsing ref introns by orthologous intron group"
	cd raw_introns_refs; \
	mkdir -p logs; \
	ls $(addprefix ../,$(INTRONREFS)) > unparsed_introns_refs.fofn; \
	/usr/bin/time -v ../../tools/parseFASTARecords.pl -f unparsed_introns_refs.fofn -m ../$(word 1,$^) -dddd 2> logs/parseFASTA.stderr > logs/parseFASTA.stdout; \
	num_spp=$$(awk 'END{print NR;}' unparsed_introns_refs.fofn); \
	for i in *_unwrapped.fasta; do num_parsed=$$(awk '/^>/{count++;}END{print count;}' $${i}); echo "$${num_parsed}"; if [[ "$${num_parsed}" -ne "$${num_spp}" ]]; then rm $${i}; fi; done | sort | uniq -c > ../$@; \
	cd ..

#Trim orthologous intron map headers and exclude requested refs:
$(SAMPLEINTRONMAP) : $(INTRONMAP)
	@echo "Trimming orthologous intron map headers for use with pseudorefs"
	../tools/trimSCOmapHeaders.awk -v "excludedcols=$(EXCLUDEDREFS)" $< > $@

#Identify orthologous introns from single copy orthogroups:
$(INTRONMAP) : $(SCOMAP) $(INTRONREFS)
	@echo "Identifying orthologous introns from single copy orthogroups"
	@echo "Using flanking exon length threshold of $(THRESHOLD)"
	ls $(INTRONREFS) > raw_introns_refs.fofn; \
	/usr/bin/time -v ../tools/findSingleCopyOrthologIntrons.pl -f raw_introns_refs.fofn -m <(awk 'BEGIN{FS="\t";OFS=FS;}NR==1{gsub("_proteome", "", $$0); print $$0;}NR>1{print $$0;}' $(word 1,$^)) -t $(THRESHOLD) -dd 2> logs/fSCOI.stderr > $@

#Trim SCO map headers and exclude requested refs:
$(SAMPLESCOMAP) : $(SCOMAP)
	@echo "Trimming single copy orthogroup map headers for use with pseudorefs"
	../tools/trimSCOmapHeaders.awk -v "excludedcols=$(EXCLUDEDREFS)" $< > $@

#Extract introns from refs:
unparsed_introns_refs/%.fasta : ../refs/%.fasta ../annotations/%.gff3 unparsed_introns_refs
	mkdir -p logs; \
	/usr/bin/time -v ../tools/extractIntronsFromGFF3.pl -d -b CDS -i $(word 1,$^) -g $(word 2,$^) -p $* > $@ 2> logs/intronome_$*.stderr

$(SUBDIRS) :
	mkdir -p $@

clean :
	for i in $(SUBDIRS); do rm -f $${i}/*.fasta $${i}/logs/* $${i}/*.fofn; [[ ! -d $${i}/logs ]] || rmdir $${i}/logs/; [[ ! -d $${i} ]] || rmdir $${i}; done
	rm -f clustalo_logs/* logs/*
	[[ ! -d clustalo_logs ]] || rmdir clustalo_logs
	[[ ! -d logs ]] || rmdir logs
	rm -f $(SAMPLESCOMAP) $(INTRONMAP) $(SAMPLEINTRONMAP) raw_introns_refs.fofn
	rm -f parse_ref_introns align_refs parse_sample_introns align_all
	rm -f samples_per_OIG_unaligned.tsv trimmed_aligned_all_number_of_seqs.tsv
