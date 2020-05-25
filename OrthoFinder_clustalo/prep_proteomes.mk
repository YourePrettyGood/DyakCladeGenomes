SHELL=/bin/bash
#Setup for OrthoFinder run
#Assumes refs and annotations are available in subdirectories
# of parent directory
#Expects the following in your PATH:
#EMBOSS transeq
#awk
#grep
#constructCDSesFromGFF3.pl
#extractSingleCopyOrthogroups.awk
#prepForOrthoFinder.awk
#remapOrthogroups.pl
#These latter three are available in the tools subdirectory of the
# Github repository.

.PHONY : SCOmap proteome clean intermediateclean

#Adjust this based on the path to OrthoFinder orthogroups output:
ORTHOFINDEROUT := OrthoFinder_proteomes/OrthoFinder/Results_Dec08/Orthogroups

REFS := $(wildcard ../refs/*.fasta)
SPECIES := $(basename $(notdir $(REFS)))
#TXOMES := $(addsuffix _txome.fna, $(SPECIES))
SUBDIRS := OrthoFinder_proteomes logs
INTERPROTEOMES := $(addsuffix _proteome.faa, $(SPECIES))
PROTEOMES := $(addprefix OrthoFinder_proteomes/, $(addsuffix _proteome.faa, $(SPECIES)))
MAPS := $(addsuffix .map, $(SPECIES))

usage :
	@echo "Usage:"
	@echo "make -f prep_proteomes.mk [task]"
	@echo "Tasks:"
	@echo "proteome -> Extracts CDSes from refs, translates, renames proteins, and generates map for renaming"
	@echo "SCOmap -> (Run after OrthoFinder) Filters Orthogroups for single-copy only, and renames genes using maps, setting up column names for use with CDS and intron alignment makefiles"
	@echo "SCOcheck -> Compares numbers of single copy orthogroups between OrthoFinder output and our awk script"
	@echo "clean -> Completely removes all created subdirectories and files"
	@echo "intermediateclean -> Removes intermediate proteome FASTAs (i.e. transeq output, not renamed)"

#OrthoFinder renamed the SingleCopyOrthologues.txt file between 2.2.7 and
# 2.3.8, so we use the 2.3.8 name here (Orthogroups_SingleCopyOrthologues.txt):
SCOcheck : $(ORTHOFINDEROUT)/Orthogroups_SingleCopyOrthologues.txt Orthogroups_SingleCopy_renamed.tsv
	printf "OrthoFinder\t" > $@; \
	awk 'END{print NR;}' $(word 1,$^) >> $@; \
	printf "awk\t" >> $@; \
	awk 'END{print NR-1;}' $(word 2,$^) >> $@; \
	printf "overlap\t" >> $@; \
	fgrep -c -f $(word 1,$^) $(word 2,$^) >> $@

SCOmap : Orthogroups_SingleCopy_renamed.tsv

Orthogroups_SingleCopy_renamed.tsv : Orthogroups_renamed.tcsv
	@echo "Filtering to retain only single-copy orthogroups"; \
	../tools/extractSingleCopyOrthogroups.awk $< > $@

#OrthoFinder renamed Orthogroups.csv to Orthogroups.tsv between 2.2.7 and
# 2.3.8, so we use the 2.3.8 version here:
Orthogroups_renamed.tcsv : combined_prots.map $(ORTHOFINDEROUT)/Orthogroups.tsv
	@echo "Renaming transcripts in Orthogroup file based on concatenated transcript renaming map"; \
	mkdir -p logs; \
	../tools/remapOrthogroups.pl -i $(word 1,$^) -g $(word 2,$^) 2> logs/remapOrthogroups.stderr > $@

combined_prots.map : $(MAPS)
	@echo "Concatenating transcript renaming maps"; \
	cat $^ > $@

proteome : $(PROTEOMES)

.SECONDARY : $(INTERPROTEOMES)

OrthoFinder_proteomes/%_proteome.faa : %_proteome.faa OrthoFinder_proteomes
	../tools/prepForOrthoFinder.awk -v "species=$*" $(word 1,$^) > $@

%_proteome.faa : %_txome.fna logs
	transeq -trim -sequence $(word 1,$^) -outseq $@ 2> logs/transeq_$*.stderr

%_txome.fna : ../refs/%.fasta ../annotations/%.gff3 logs
	../tools/constructCDSesFromGFF3.pl -l -i $(word 1,$^) -g $(word 2,$^) -p $* > $@ 2> logs/txome_$*.stderr

$(SUBDIRS) :
	mkdir -p $@

clean :
	rm -f $(PROTEOMES) $(INTERPROTEOMES) logs/*.stderr *.map
	rm -f Orthogroups_renamed.tcsv Orthogroups_SingleCopy_renamed.tsv SCOcheck
	[[ ! -d OrthoFinder_proteomes ]] || rmdir OrthoFinder_proteomes
	[[ ! -d logs ]] || rmdir logs

intermediateclean :
	rm -f $(INTERPROTEOMES)
