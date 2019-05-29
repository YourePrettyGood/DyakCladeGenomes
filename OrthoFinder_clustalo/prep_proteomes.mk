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
#remapOrthogroups.pl
#These latter three are available in the tools subdirectory of the
# Github repository.

.PHONY : SCOmap proteome clean intermediateclean

#Adjust this based on the path to OrthoFinder output:
ORTHOFINDEROUT := OrthoFinder_proteomes/Results_May16

REFS := $(wildcard ../refs/*.fasta)
SPECIES := $(basename $(notdir $(REFS)))
#SPECIES := DmelISO1 DsanSTOCAGO1482 DteiGT53w DyakNY73PB DyakTai18E2
#TXOMES := $(addsuffix _txome.fna, $(SPECIES))
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

SCOcheck : $(ORTHOFINDEROUT)/SingleCopyOrthogroups.txt Orthogroups_SingleCopy_renamed.tsv
	printf "OrthoFinder\t" > $@; \
	awk 'END{print NR;}' $(word 1,$^) >> $@; \
	printf "awk\t" >> $@; \
	awk 'END{print NR-1;}' $(word 2,$^) >> $@; \
	printf "overlap\t" >> $@; \
	fgrep -c -f $(word 1,$^) $(word 2,$^) >> $@

SCOmap : Orthogroups_SingleCopy_renamed.tsv

Orthogroups_SingleCopy_renamed.tsv : Orthogroups_renamed.tcsv
	@echo "Filtering to retain only single-copy orthogroups"; \
	extractSingleCopyOrthogroups.awk $< > $@

Orthogroups_renamed.tcsv : combined_prots.map $(ORTHOFINDEROUT)/Orthogroups.csv
	@echo "Renaming transcripts in Orthogroup file based on concatenated transcript renaming map"; \
	mkdir -p logs; \
	remapOrthogroups.pl -i $(word 1,$^) -g $(word 2,$^) 2> logs/remapOrthogroups.stderr > $@

combined_prots.map : $(MAPS)
	@echo "Concatenating transcript renaming maps"; \
	cat $^ > $@

proteome : $(PROTEOMES)

.SECONDARY : $(INTERPROTEOMES)

OrthoFinder_proteomes/%_proteome.faa : %_proteome.faa
	awk -v "species=$*" 'BEGIN{i=0; mapfn=species".map"; speciespattern=">"species"_";}!/^>/{print;}/^>/{print ">"species"_"i; txid=$$0; sub(speciespattern, "", txid); sub(/_1$$/, "", txid); print txid"\t"species"_"i > mapfn; i+=1;}' $< > $@

%_proteome.faa : %_txome.fna
	mkdir -p OrthoFinder_proteomes
	transeq -trim -sequence $< -outseq $@ 2> logs/transeq_$*.stderr

%_txome.fna : ../refs/%.fasta ../annotations/%.gff3
	mkdir -p logs
	constructCDSesFromGFF3.pl -l -i $(word 1,$^) -g $(word 2,$^) -p $* > $@ 2> logs/txome_$*.stderr

clean :
	rm -f $(PROTEOMES) $(INTERPROTEOMES) logs/*.stderr *.map
	rm -f Orthogroups_renamed.tcsv Orthogroups_SingleCopy_renamed.tsv SCOcheck
	[[ ! -d OrthoFinder_proteomes ]] || rmdir OrthoFinder_proteomes
	[[ ! -d logs ]] || rmdir logs

intermediateclean :
	rm -f $(INTERPROTEOMES)
