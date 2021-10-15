SHELL=/bin/bash

#Calculate genome-wide pi (at all sites and 4-fold sites only)
#Also generate windowed estimates of pi, and depth
#Assumes refs and annotations are available in subdirectories
# of parent directory, plus pseudoreferences
#For now, we separately extracted per-site depth from VCFs
# and placed them in the depth subdirectory

#Expects the following in your PATH:
#bedtools
#fasta_formatter
#Expects these Awk, C++, and Perl scripts in the tools subdirectory:
#calculatePolymorphism
#CDStoGenomicIntervals.pl
#codingSitesByDegeneracy.pl
#constructCDSesFromGFF3.pl
#decompressStats.pl
#nonOverlappingWindows
#subsetVCFstats.pl
#extractOGlocation.awk
#selectMajorArms.awk
#OGIDsFromBED.awk
#These Awk, C++, and Perl scripts are available in the tools subdirectory
# of the Github repository.
#Also expects the following R script in the tools subdirectory:
#IdentifyLowPiRegions.R

#Options to adjust:
#Window size (bp):
WINDOW := 100000
#Ref(s) to exclude from analysis:
EXCLUDEDREFS := Dyak_Tai18E2
#Treat pseudorefs as inbred, so only use one allele from each:
#Leave empty if all are outbred diploid, use -i if any are inbred
INBRED := -i
#PRNG seed to use for selecting random allele if inbred:
#Do not leave blank
SEED := 42
#Infimum non-N fraction for a site to be included in the windowed average:
INFNONN := 0.5
#Chromosome arms to check for low-pi regions:
ARMS := X 2L 2R 3L 3R

#DO NOT CHANGE THE BELOW:
#Subdirectories that need to be created:
SUBDIRS := CDSes BEDs logs
#List of pseudorefs available:
PSEUDOREFS := $(wildcard ../pseudorefs/*.fasta)
#List of samples used:
SAMPLES := $(patsubst %.fasta,%,$(notdir $(PSEUDOREFS)))

#List of SCOs in normal-pi regions in all three members of the Dyak clade:
CLADEINTERSECTION := Dsan_Dtei_Dyak_normal_pi_SCOs_intersection.txt
#References inferred from directory above the current:
REFS := $(filter-out $(EXCLUDEDREFS),$(basename $(notdir $(wildcard ../refs/*.fasta))))
#Per-reference orthogroup location files:
OGSPPLOCS := $(addsuffix _OG_locs.tsv,$(REFS))
#Per-reference orthogroup location BEDs:
OGLOCBEDS := $(addprefix BEDs/,$(addsuffix _SCOs.bed,$(REFS)))
#Line-specific BEDs of 4-fold sites:
REFBEDS := $(addprefix BEDs/,$(addsuffix _line_4fold_genomic.bed,$(REFS)))
#Species that the references represent:
SPECIES := $(foreach ref,$(REFS),$(word 1,$(subst _, ,$(ref))))
#Species-specific FAIs:
REFFAIS := $(addprefix FAIs/,$(addsuffix .fasta.fai,$(SPECIES)))
#Genome-wide per-site pi estimate files:
GWPI := $(addsuffix _genomewide_pi.tsv.gz,$(SPECIES))
#Windowed genome-wide pi estimate files:
WINDOWEDGWPI := $(addsuffix _genomewide_pi_w$(WINDOW).tsv.gz,$(SPECIES))
#4-fold per-site pi estimate files:
FFPI := $(addsuffix _4fold_pi.tsv.gz,$(SPECIES))
#Windowed 4-fold pi estimate files:
WINDOWEDFFPI := $(addsuffix _4fold_pi_w$(WINDOW).tsv.gz,$(SPECIES))
#4-fold depth files:
FFDP := $(addprefix depth/,$(addsuffix _4fold_depth.tsv.gz,$(SAMPLES)))
#Windowed 4-fold depth files:
WINDOWEDFFDP := $(addprefix depth/,$(addsuffix _4fold_depth_w$(WINDOW).tsv.gz,$(SAMPLES)))
#Low-pi BED files:
LOWPIBED := $(addprefix BEDs/,$(addsuffix _low_pi_regions.bed,$(SPECIES)))
#Normal-pi BED files:
NORMALPIBED := $(addprefix BEDs/,$(addsuffix _normal_pi_regions.bed,$(SPECIES)))
#Per-reference orthogroup location BEDs with only species as prefix:
OGLOCBEDSSPP := $(addprefix BEDs/,$(addsuffix _SCOs.bed,$(SPECIES)))
#Normal-pi SCO BED files:
NORMALPISCO := $(addprefix BEDs/,$(addsuffix _normal_pi_SCOs.bed,$(SPECIES)))
#Normal-pi SCO OGID lists:
NORMALPISCOTXT := $(addprefix BEDs/,$(addsuffix _normal_pi_SCOs.txt,$(SPECIES)))

.PHONY : all orthogroup_locations low_pi clean usage

.SECONDARY : $(LOWPIBED) $(NORMALPIBED) $(NORMALPISCO) $(NORMALPISCOTXT) $(OGLOCBEDS) $(GWPI) $(WINDOWEDGWPI) $(FFPI) $(WINDOWEDFFPI) $(REFFAIS) $(REFBEDS)

usage :
	@echo "Usage:"
	@echo "make -f calculate_pi.mk [task]"
	@echo "Tasks:"
	@echo "all -> Calculate per-site and windowed genome-wide and 4-fold-only pi for each species"
	@echo "orthogroup_locations -> Find SCO locations"
	@echo "low_pi -> Identify low-pi regions and partition SCOs by this"
	@echo "clean -> Clean up all output files and directories"

#all : $(GWPI) $(WINDOWEDGWPI) $(FFPI) $(WINDOWEDFFPI) $(WINDOWEDFFDP)
all : $(GWPI) $(WINDOWEDGWPI) $(FFPI) $(WINDOWEDFFPI)

low_pi : $(CLADEINTERSECTION) $(NORMALPISCOTXT)

#This rule is custom for the Dyak clade, so it will need to be edited
# to be general:
$(CLADEINTERSECTION) : BEDs/Dsan_normal_pi_SCOs.txt BEDs/Dtei_normal_pi_SCOs.txt BEDs/Dyak_normal_pi_SCOs.txt
	@echo "Taking intersection of Dyak clade normal-pi OGID lists"
	sort $^ | uniq -c | awk '$$1==3{print $$2;}' > $@

$(NORMALPISCOTXT) : BEDs/%_normal_pi_SCOs.txt : BEDs/%_normal_pi_SCOs.bed
	@echo "Making OGID lists for normal-pi SCOs in $*"
	../tools/OGIDsFromBED.awk $< > $@

$(NORMALPISCO) : BEDs/%_normal_pi_SCOs.bed : BEDs/%_SCOs.bed BEDs/%_normal_pi_regions.bed
	@echo "Extracting SCOs with normal pi in $*"
	../tools/subsetVCFstats.pl -i $(word 1,$^) -b $(word 2,$^) > $@

#As a mediocre solution, let's force all of the per-ref BEDs to exist,
# but only link one of them:
$(OGLOCBEDSSPP) : BEDs/%_SCOs.bed : $(OGLOCBEDS)
	cd BEDs/; \
	ln -sf $*_*_SCOs.bed `basename $@`; \
	cd ../

$(NORMALPIBED) : BEDs/%_normal_pi_regions.bed : BEDs/%_low_pi_regions.bed FAIs/%.genome
	@echo "Selecting normal-pi regions for major chromosome arms of $*"
	bedtools complement -i <(sort -k1,1 -k2,2n $(word 1,$^)) -g <(sort -k1,1 $(word 2,$^)) | \
	../tools/selectMajorArms.awk > $@

$(LOWPIBED) : BEDs/%_low_pi_regions.bed : FAIs/%.fasta.fai %_genomewide_pi_w$(WINDOW).tsv.gz logs BEDs
	@echo "Identifying low-pi regions for major chromosome arms of $*"
	../tools/IdentifyLowPiRegions.R $(word 1,$^) $(word 2,$^) $* $(ARMS) 2> $(word 3,$^)/ILPR_$*.stderr > $@

FAIs/%.genome : FAIs/%.fasta.fai
	cut -f1,2 $< > $@

#Windowed depth calculations:
#$(WINDOWEDFFDP) : depth/%_4fold_depth_w$(WINDOW).tsv.gz : depth/%_4fold_depth.tsv.gz
#	@echo "Calculating windowed depth for $*"
#	../tools/nonOverlappingWindows -w $(WINDOW) -f $(INFNONN) -u -i <(gzip -dc $<) 2> logs/nOW_$*_4fold_depth_w$(WINDOW).stderr | gzip -9 > $@

#4-fold depth subsetting:
#$(FFDP) : %_4fold_depth.tsv.gz : %_DP_persite.tsv.gz BEDs/%_spp_4fold_genomic.bed FAIs/%.fasta.fai
#	@echo "Subsetting 4-fold sites from genome-wide depth for $*"
#	mkdir -p logs; \
#	../tools/subsetVCFstats.pl -i <(gzip -dc $(word 1,$^)) -b $(word 2,$^) 2> logs/sVs_$*_4fold_depth.stderr | ../tools/decompressStats.pl -f $(word 3,$^) -u 2> logs/dS_$*_4fold_depth.stderr | gzip -9 > $@

#Windowed pi calculations:
%_pi_w$(WINDOW).tsv.gz : %_pi.tsv.gz
	@echo "Calculating windowed pi for $*"
	mkdir -p logs; \
	../tools/nonOverlappingWindows -w $(WINDOW) -f $(INFNONN) -u -i <(gzip -dc $<) 2> logs/nOW_$*_w$(WINDOW).stderr | gzip -9 > $@

#4-fold genome-wide pi subsetting:
%_4fold_pi.tsv.gz : %_genomewide_pi.tsv.gz BEDs/%_spp_4fold_genomic.bed FAIs/%.fasta.fai
	@echo "Subsetting 4-fold sites from genome-wide pi for $*"
	mkdir -p logs; \
	../tools/subsetVCFstats.pl -i <(gzip -dc $(word 1,$^)) -b $(word 2,$^) 2> logs/sVs_$*_4fold.stderr | ../tools/decompressStats.pl -f $(word 3,$^) -u 2> logs/dS_$*_4fold.stderr | gzip -9 > $@

$(REFFAIS) : FAIs/%.fasta.fai : FAIs/%.fasta
	@echo "Indexing reference $* with samtools"
	samtools faidx $<

FAIs/%.fasta : FAIs
	cd FAIs/; \
	ln -sf $*_*.fasta $*.fasta; \
	cd ..

FAIs : 
	mkdir -p $@; \
	ln -sf $(addprefix ../../refs/,$(addsuffix .fasta,$(REFS))) FAIs/

#Genome-wide pi calculations:
%_genomewide_pi.tsv.gz : %_pseudorefs.fofn
	@echo "Calculating per-base pi genomewide for $*"
	mkdir -p logs; \
	../tools/calculatePolymorphism $(INBRED) -u -p $(SEED) -f $< 2> logs/cP_$*_genomewide.stderr | gzip -9 > $@

%_pseudorefs.fofn : $(PSEUDOREFS)
	ls ../pseudorefs/$*_*.fasta > $@

#Make symlinks from ref line-specific BEDs to species BEDs:
BEDs/%_spp_4fold_genomic.bed : $(REFBEDS)
	BED=`basename $@`; \
	cd BEDs; \
	ln -sf $*_*_line_4fold_genomic.bed $${BED}; \
	cd ..

#Define 4-fold sites in genome space:
BEDs/%_line_4fold_genomic.bed : BEDs/%_4fold_CDS.bed ../annotations/%.gff3
	@echo "Converting $* 4-fold site coordinates into genome-space"
	mkdir -p logs; \
	../tools/CDStoGenomicIntervals.pl -i $(word 1,$^) -g $(word 2,$^) 2> logs/CtGI_$*_4fold.stderr | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > $@

#Define 4-fold sites in CDS space:
BEDs/%_4fold_CDS.bed : CDSes/%_CDSes.fasta BEDs
	@echo "Identifying $* 4-fold sites in CDS coordinates"
	mkdir -p logs; \
	../tools/codingSitesByDegeneracy.pl -f 4 -i $(word 1,$^) 2> logs/cSBD_$*_4fold.stderr | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > $@

#Make CDSes needed to map out degeneracy of sites:
CDSes/%_CDSes.fasta : ../refs/%.fasta ../annotations/%.gff3 CDSes
	@echo "Extracting CDSes for $*"
	mkdir -p logs; \
	../tools/constructCDSesFromGFF3.pl -i $(word 1,$^) -g $(word 2,$^) -l 2> logs/cCFG_$*.stderr > $@

orthogroup_locations : $(OGSPPLOCS)
	cat $^ | sort -k5,5V > $@

$(OGLOCBEDS) : BEDs/%_SCOs.bed : %_OG_locs.tsv BEDs
	@echo "Making BED of SCO locations for $*"
	../tools/SCOlocationBED.awk -v "spid=$*_" $(word 1,$^) | \
	sort -k1,1 -k2,2n -k3,3n > $@

$(OGSPPLOCS) : %_OG_locs.tsv : ../OrthoFinder_clustalo/Orthogroups_SingleCopy_renamed.tsv ../annotations/%.gff3
	@echo "Identifying $* locations for single-copy orthogroups (SCOs)"
	../tools/extractOGlocation.awk -v "refid=$*_proteome" $(word 1,$^) $(word 2,$^) > $@

$(SUBDIRS) :
	mkdir -p $@

clean :
	for i in $(SUBDIRS); do rm -f $${i}/*.fasta $${i}/*.fai $${i}/*.bed $${i}/*.txt $${i}/*.stderr; [[ ! -d $${i} ]] || rmdir $${i}; done
	rm -f $(GWPI) $(WINDOWEDGWPI) $(FFPI) $(WINDOWEDFFPI) $(OGSPPLOCS)
	rm -rf FAIs/
	rm -f *.fofn
