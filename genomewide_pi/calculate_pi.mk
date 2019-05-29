SHELL=/bin/bash

#Calculate genome-wide pi (at all sites and 4-fold sites only)
#Also generate windowed estimates of pi
#Assumes refs and annotations are available in subdirectories
# of parent directory, plus pseudoreferences

#Expects the following in your PATH:
#bedtools
#fasta_formatter
#calculatePolymorphism
#CDStoGenomicIntervals.pl
#codingSitesByDegeneracy.pl
#constructCDSesFromGFF3.pl
#decompressStats.pl
#nonOverlappingWindows
#subsetVCFstats.pl
#These Perl and C++ scripts are available in the tools subdirectory
# of the Github repository.

#Options to adjust:
#Window size (bp):
WINDOW := 100000
#Ref(s) to exclude from analysis:
EXCLUDEDREFS := Dyak_Tai18E2 Dmel_ISO1
#Treat pseudorefs as inbred, so only use one allele from each:
#Leave empty if all are outbred diploid, use -i if any are inbred
INBRED := -i
#PRNG seed to use for selecting random allele if inbred:
#Do not leave blank
SEED := 42

#DO NOT CHANGE THE BELOW:
#Subdirectories that need to be created:
SUBDIRS := CDSes BEDs logs
#List of pseudorefs available:
PSEUDOREFS := $(wildcard ../pseudorefs/*.fasta)

#References inferred from directory above the current:
REFS := $(filter-out $(EXCLUDEDREFS),$(basename $(notdir $(wildcard ../refs/*.fasta))))
#Line-specific BEDs of 4-fold sites:
REFBEDS := $(addprefix BEDs/,$(addsuffix _line_4fold_genomic.bed,$(REFS)))
#Species that the references represent:
SPECIES := $(foreach ref,$(REFS),$(word 1,$(subst _, ,$(ref))))
#Genome-wide per-site pi estimate files:
GWPI := $(addsuffix _genomewide_pi.tsv.gz,$(SPECIES))
#Windowed genome-wide pi estimate files:
WINDOWEDGWPI := $(addsuffix _genomewide_pi_w$(WINDOW).tsv.gz,$(SPECIES))
#4-fold per-site pi estimate files:
FFPI := $(addsuffix _4fold_pi.tsv.gz,$(SPECIES))
#Windowed 4-fold pi estimate files:
WINDOWEDFFPI := $(addsuffix _4fold_pi_w$(WINDOW).tsv.gz,$(SPECIES))


.PHONY : all clean usage

.SECONDARY : $(GWPI) $(WINDOWEDGWPI) $(FFPI) $(WINDOWEDFFPI)

usage :
	@echo "Usage:"
	@echo "make -f calculate_pi.mk [task]"
	@echo "Tasks:"
	@echo "all -> Calculate per-site and windowed genome-wide and 4-fold-only pi for each species"
	@echo "clean -> Clean up all output files and directories"

all : $(GWPI) $(WINDOWEDGWPI) $(FFPI) $(WINDOWEDFFPI)

#Windowed pi calculations:
%_pi_w$(WINDOW).tsv.gz : %_pi.tsv.gz
	@echo "Calculating windowed pi for $*"
	mkdir -p logs; \
	nonOverlappingWindows -w $(WINDOW) -a -u -i <(zcat $<) 2> logs/nOW_$*_w$(WINDOW).stderr | gzip -9 > $@

#4-fold genome-wide pi subsetting:
%_4fold_pi.tsv.gz : %_genomewide_pi.tsv.gz BEDs/%_spp_4fold_genomic.bed
	@echo "Subsetting 4-fold sites from genome-wide pi for $*"
	mkdir -p logs; \
	subsetVCFstats.pl -i <(gzip -dc $(word 1,$^)) -b $(word 2,$^) 2> logs/sVs_$*_4fold.stderr | gzip -9 > $@

#Genome-wide pi calculations:
%_genomewide_pi.tsv.gz : %_pseudorefs.fofn
	@echo "Calculating genome-wide pi for $*"
	mkdir -p logs; \
	calculatePolymorphism $(INBRED) -u -p $(SEED) -f $< 2> logs/cP_$*_genomewide.stderr | gzip -9 > $@

%_pseudorefs.fofn : $(PSEUDOREFS)
	ls ../pseudorefs/$*_*.fasta > $@

#Make symlinks from ref line-specific BEDs to species BEDs:
BEDs/%_spp_4fold_genomic.bed : $(REFBEDS)
	BED=`basename $@`; \
	cd BEDs; \
	ln -s $*_*_line_4fold_genomic.bed $${BED}; \
	cd ..

#Define 4-fold sites in genome space:
BEDs/%_line_4fold_genomic.bed : BEDs/%_4fold_CDS.bed ../annotations/%.gff3
	@echo "Converting $* 4-fold site coordinates into genome-space"
	mkdir -p logs; \
	CDStoGenomicIntervals.pl -i $(word 1,$^) -g $(word 2,$^) 2> logs/CtGI_$*_4fold.stderr | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > $@

#Define 4-fold sites in CDS space:
BEDs/%_4fold_CDS.bed : CDSes/%_CDSes.fasta BEDs
	@echo "Identifying $* 4-fold sites in CDS coordinates"
	mkdir -p logs; \
	codingSitesByDegeneracy.pl -f 4 -i $(word 1,$^) 2> logs/cSBD_$*_4fold.stderr | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > $@

#Make CDSes needed to map out degeneracy of sites:
CDSes/%_CDSes.fasta : ../refs/%.fasta ../annotations/%.gff3 CDSes
	@echo "Extracting CDSes for $*"
	mkdir -p logs; \
	constructCDSesFromGFF3.pl -i $(word 1,$^) -g $(word 2,$^) -l 2> logs/cCFG_$*.stderr > $@

$(SUBDIRS) :
	mkdir -p $@

clean :
	for i in $(SUBDIRS); do rm -f $${i}/*.fasta $${i}/*.bed $${i}/*.stderr; [[ ! -d $${i} ]] || rmdir $${i}; done
	rm -f $(GWPI) $(WINDOWEDGWPI) $(FFPI) $(WINDOWEDFFPI)
	rm -f *.fofn