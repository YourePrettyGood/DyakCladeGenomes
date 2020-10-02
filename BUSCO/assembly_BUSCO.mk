SHELL=/bin/bash

#Runs BUSCO v3 on the assemblies, and generates a BUSCO summary plot
#Assumes refs are available in subdirectory of parent directory
#For transcriptome and/or proteome target, assumes annotations
# are available in subdirectory of parent directory

#Expects the following in your PATH:
#run_BUSCO.py
#generate_plot.py
#All BUSCO v3 dependencies (i.e. NCBI BLAST+, HMMer 3.1b2, Augustus 3.x)
#generate_plot.py R dependencies (i.e. ggplot2)
#For transcriptome and proteome BUSCO, requires the following in your PATH:
#constructCDSesFromGFF3.pl
#EMBOSS transeq
#The Perl script can be found in the tools directory of the Github repo

#Options to adjust:
LINEAGE := diptera
LINEAGEPATH := /home/pfreilly/Downloads/Bioinformatics/Assembly/busco/diptera_odb9/
#Old refs that don't have annotations (space-separated list of prefixes):
NOANNOT := Dyak_Miller

#DO NOT CHANGE THE BELOW:
#Number of cores to use per BUSCO job (use 8, as kfold for
# optimize_augustus.pl is set to 8, so anything different will be
# inefficient and produce warnings):
NUMCORES := 8
#Subdirectories to create:
SUBDIRS := logs BUSCO_summaries BUSCO_txome_summaries BUSCO_proteome_summaries
#References inferred from directory above the current:
REFS := $(basename $(notdir $(wildcard ../refs/*.fasta)))
#Old reference(s) (not PacBio one(s)), inferred from "old" subdirectory of
# "refs" directory above current:
OLDREFS := $(basename $(notdir $(wildcard ../refs/old/*.fasta)))
#Unpublished reference, i.e. Dsim w501 v3 unpublished:
UNPUBREFS := $(basename $(notdir $(wildcard ../refs/Dsim_w501_v3_unpublished/*.fasta)))
#BUSCO result files:
BUSCOS := $(addsuffix _buscos.tsv,$(REFS))
#BUSCO result file(s) for old references (not our PacBio ones):
EXTRABUSCOS := $(addsuffix _buscos.tsv,$(OLDREFS))
#BUSCO result files for unpublished reference:
UNPUBBUSCOS := $(addsuffix _unpub_buscos.tsv,$(UNPUBREFS))
#CDS FASTAs:
CDSES := $(addsuffix _CDSes.fasta,$(REFS))
#CDS FASTA(s) for old references (not our PacBio ones):
EXTRACDSES := $(addsuffix _CDSes.fasta,$(OLDREFS))
#CDS FASTA for unpublished reference:
UNPUBCDSES := $(addsuffix _unpub_CDSes.fasta,$(UNPUBREFS))
#Transcriptome BUSCO result files:
TRANBUSCOS := $(addsuffix _txome_buscos.tsv,$(REFS))
#Transcriptome BUSCO result file(s) for old references (not our PacBio ones):
EXTRATRANBUSCOS := $(addsuffix _txome_buscos.tsv,$(filter-out $(NOANNOT),$(OLDREFS)))
#Transcriptome BUSCO result files for unpublished reference:
UNPUBTRANBUSCOS := $(addsuffix _unpub_txome_buscos.tsv,$(filter-out $(NOANNOT),$(UNPUBREFS)))
#Proteome BUSCO result files:
PROTBUSCOS := $(addsuffix _proteome_buscos.tsv,$(REFS))
#Proteome BUSCO result file(s) for old references (not our PacBio ones):
EXTRAPROTBUSCOS := $(addsuffix _proteome_buscos.tsv,$(filter-out $(NOANNOT),$(OLDREFS)))
#Proteome BUSCO result files for unpublished reference:
UNPUBPROTBUSCOS := $(addsuffix _unpub_proteome_buscos.tsv,$(filter-out $(NOANNOT),$(UNPUBREFS)))

.PHONY : genome_plot transcriptome_plot proteome_plot run_unpublished clean usage

.SECONDARY : $(BUSCOS) $(EXTRABUSCOS) genome transcriptome proteome run_unpublished

usage :
	@echo "Usage:"
	@echo "make -f assembly_BUSCO.mk [task]"
	@echo "Tasks:"
	@echo "genome_plot -> Make a combined BUSCO plot from all genomes"
	@echo "genome -> Run BUSCO genome analysis with 8 cores per job"
	@echo "transcriptome_plot -> Make a combined BUSCO plot from all transcriptomes"
	@echo "transcriptome -> Run BUSCO transcriptome analysis"
	@echo "proteome_plot -> Make a combined BUSCO plot from all proteomes"
	@echo "proteome -> Run BUSCO proteome analysis"
	@echo "run_unpublished -> Run all BUSCO analyses for unpublished Dsim w501 reference"
	@echo "clean -> Clean up all output files and directories"

run_unpublished : $(UNPUBBUSCOS) $(UNPUBTRANBUSCOS) $(UNPUBPROTBUSCOS)

genome_plot : genome
	generate_plot.py -wd BUSCO_summaries

genome : $(BUSCOS) $(EXTRABUSCOS)
	cat $^ > $@

#Run BUSCO v3 on each assembly, using 8 cores per job:
$(UNPUBBUSCOS) : %_unpub_buscos.tsv : ../refs/Dsim_w501_v3_unpublished/%.fasta logs BUSCO_summaries
	@echo "Running BUSCO v3 on $*_unpub"
	/usr/bin/time -v run_BUSCO.py -i $(word 1,$^) -l $(LINEAGEPATH) -m genome -c $(NUMCORES) --long -z -o $*_unpub_BUSCOv3_genome_$(LINEAGE) 2> $(word 2,$^)/$*_unpub_BUSCOv3_genome_$(LINEAGE).stderr > $(word 2,$^)/$*_unpub_BUSCOv3_genome_$(LINEAGE).stdout; \
	tail -n10 $(word 2,$^)/$*_unpub_BUSCOv3_genome_$(LINEAGE).stdout | head -n7 | awk 'BEGIN{printf "%s", "$*";}NR>1{printf "\t%s", $$2;}END{print "";}' > $@; \
	cp run_$*_unpub_BUSCOv3_genome_$(LINEAGE)/short_summary_$*_unpub_BUSCOv3_genome_$(LINEAGE).txt $(word 3,$^)/short_summary_$*_unpub.txt

$(EXTRABUSCOS) : %_buscos.tsv : ../refs/old/%.fasta logs BUSCO_summaries
	@echo "Running BUSCO v3 on $*"
	/usr/bin/time -v run_BUSCO.py -i $(word 1,$^) -l $(LINEAGEPATH) -m genome -c $(NUMCORES) --long -z -o $*_BUSCOv3_genome_$(LINEAGE) 2> $(word 2,$^)/$*_BUSCOv3_genome_$(LINEAGE).stderr > $(word 2,$^)/$*_BUSCOv3_genome_$(LINEAGE).stdout; \
	tail -n10 $(word 2,$^)/$*_BUSCOv3_genome_$(LINEAGE).stdout | head -n7 | awk 'BEGIN{printf "%s", "$*";}NR>1{printf "\t%s", $$2;}END{print "";}' > $@; \
	cp run_$*_BUSCOv3_genome_$(LINEAGE)/short_summary_$*_BUSCOv3_genome_$(LINEAGE).txt $(word 3,$^)/short_summary_$*.txt

$(BUSCOS) : %_buscos.tsv : ../refs/%.fasta logs BUSCO_summaries
	@echo "Running BUSCO v3 on $*"
	/usr/bin/time -v run_BUSCO.py -i $(word 1,$^) -l $(LINEAGEPATH) -m genome -c $(NUMCORES) --long -z -o $*_BUSCOv3_genome_$(LINEAGE) 2> $(word 2,$^)/$*_BUSCOv3_genome_$(LINEAGE).stderr > $(word 2,$^)/$*_BUSCOv3_genome_$(LINEAGE).stdout; \
	tail -n10 $(word 2,$^)/$*_BUSCOv3_genome_$(LINEAGE).stdout | head -n7 | awk 'BEGIN{printf "%s", "$*";}NR>1{printf "\t%s", $$2;}END{print "";}' > $@; \
	cp run_$*_BUSCOv3_genome_$(LINEAGE)/short_summary_$*_BUSCOv3_genome_$(LINEAGE).txt $(word 3,$^)/short_summary_$*.txt

transcriptome_plot : transcriptome
	generate_plot.py -wd BUSCO_txome_summaries

transcriptome : $(TRANBUSCOS) $(EXTRATRANBUSCOS)
	cat $^ > $@

#Run BUSCO v3 on each transcriptome, using 8 cores per job:
$(UNPUBTRANBUSCOS) : %_unpub_txome_buscos.tsv : %_CDSes.fasta logs BUSCO_txome_summaries
	@echo "Running BUSCO v3 transcriptome on $*_unpub"
	/usr/bin/time -v run_BUSCO.py -i $(word 1,$^) -l $(LINEAGEPATH) -m transcriptome -c $(NUMCORES) -z -o $*_unpub_BUSCOv3_txome_$(LINEAGE) 2> $(word 2,$^)/$*_unpub_BUSCOv3_txome_$(LINEAGE).stderr > $(word 2,$^)/$*_unpub_BUSCOv3_txome_$(LINEAGE).stdout; \
	tail -n10 $(word 2,$^)/$*_unpub_BUSCOv3_txome_$(LINEAGE).stdout | head -n7 | awk 'BEGIN{printf "%s", "$*";}NR>1{printf "\t%s", $$2;}END{print "";}' > $@; \
	cp run_$*_unpub_BUSCOv3_txome_$(LINEAGE)/short_summary_$*_unpub_BUSCOv3_txome_$(LINEAGE).txt $(word 3,$^)/short_summary_$*_unpub.txt

$(EXTRATRANBUSCOS) : %_txome_buscos.tsv : %_CDSes.fasta logs BUSCO_txome_summaries
	@echo "Running BUSCO v3 transcriptome on $*"
	/usr/bin/time -v run_BUSCO.py -i $(word 1,$^) -l $(LINEAGEPATH) -m transcriptome -c $(NUMCORES) -z -o $*_BUSCOv3_txome_$(LINEAGE) 2> $(word 2,$^)/$*_BUSCOv3_txome_$(LINEAGE).stderr > $(word 2,$^)/$*_BUSCOv3_txome_$(LINEAGE).stdout; \
	tail -n10 $(word 2,$^)/$*_BUSCOv3_txome_$(LINEAGE).stdout | head -n7 | awk 'BEGIN{printf "%s", "$*";}NR>1{printf "\t%s", $$2;}END{print "";}' > $@; \
	cp run_$*_BUSCOv3_txome_$(LINEAGE)/short_summary_$*_BUSCOv3_txome_$(LINEAGE).txt $(word 3,$^)/short_summary_$*.txt

$(TRANBUSCOS) : %_txome_buscos.tsv : %_CDSes.fasta logs BUSCO_txome_summaries
	@echo "Running BUSCO v3 transcriptome on $*"
	/usr/bin/time -v run_BUSCO.py -i $(word 1,$^) -l $(LINEAGEPATH) -m transcriptome -c $(NUMCORES) -z -o $*_BUSCOv3_txome_$(LINEAGE) 2> $(word 2,$^)/$*_BUSCOv3_txome_$(LINEAGE).stderr > $(word 2,$^)/$*_BUSCOv3_txome_$(LINEAGE).stdout; \
	tail -n10 $(word 2,$^)/$*_BUSCOv3_txome_$(LINEAGE).stdout | head -n7 | awk 'BEGIN{printf "%s", "$*";}NR>1{printf "\t%s", $$2;}END{print "";}' > $@; \
	cp run_$*_BUSCOv3_txome_$(LINEAGE)/short_summary_$*_BUSCOv3_txome_$(LINEAGE).txt $(word 3,$^)/short_summary_$*.txt

proteome_plot : proteome
	generate_plot.py -wd BUSCO_proteome_summaries

proteome : $(PROTBUSCOS) $(EXTRAPROTBUSCOS)
	cat $^ > $@

#Run BUSCO v3 on each proteome, using 8 cores per job:
$(UNPUBPROTBUSCOS) : %_unpub_proteome_buscos.tsv : %_unpub_proteome.fasta logs BUSCO_proteome_summaries
	@echo "Running BUSCO v3 proteome on $*_unpub"
	/usr/bin/time -v run_BUSCO.py -i $(word 1,$^) -l $(LINEAGEPATH) -m proteins -c $(NUMCORES) -z -o $*_unpub_BUSCOv3_proteome_$(LINEAGE) 2> $(word 2,$^)/$*_unpub_BUSCOv3_proteome_$(LINEAGE).stderr > $(word 2,$^)/$*_unpub_BUSCOv3_proteome_$(LINEAGE).stdout; \
	tail -n10 $(word 2,$^)/$*_unpub_BUSCOv3_proteome_$(LINEAGE).stdout | head -n7 | awk 'BEGIN{printf "%s", "$*";}NR>1{printf "\t%s", $$2;}END{print "";}' > $@; \
	cp run_$*_unpub_BUSCOv3_proteome_$(LINEAGE)/short_summary_$*_unpub_BUSCOv3_proteome_$(LINEAGE).txt $(word 3,$^)/short_summary_$*_unpub.txt

$(EXTRAPROTBUSCOS) : %_proteome_buscos.tsv : %_proteome.fasta logs BUSCO_proteome_summaries
	@echo "Running BUSCO v3 proteome on $*"
	/usr/bin/time -v run_BUSCO.py -i $(word 1,$^) -l $(LINEAGEPATH) -m proteins -c $(NUMCORES) -z -o $*_BUSCOv3_proteome_$(LINEAGE) 2> $(word 2,$^)/$*_BUSCOv3_proteome_$(LINEAGE).stderr > $(word 2,$^)/$*_BUSCOv3_proteome_$(LINEAGE).stdout; \
	tail -n10 $(word 2,$^)/$*_BUSCOv3_proteome_$(LINEAGE).stdout | head -n7 | awk 'BEGIN{printf "%s", "$*";}NR>1{printf "\t%s", $$2;}END{print "";}' > $@; \
	cp run_$*_BUSCOv3_proteome_$(LINEAGE)/short_summary_$*_BUSCOv3_proteome_$(LINEAGE).txt $(word 3,$^)/short_summary_$*.txt

$(PROTBUSCOS) : %_proteome_buscos.tsv : %_proteome.fasta logs BUSCO_proteome_summaries
	@echo "Running BUSCO v3 proteome on $*"
	/usr/bin/time -v run_BUSCO.py -i $(word 1,$^) -l $(LINEAGEPATH) -m proteins -c $(NUMCORES) -z -o $*_BUSCOv3_proteome_$(LINEAGE) 2> $(word 2,$^)/$*_BUSCOv3_proteome_$(LINEAGE).stderr > $(word 2,$^)/$*_BUSCOv3_proteome_$(LINEAGE).stdout; \
	tail -n10 $(word 2,$^)/$*_BUSCOv3_proteome_$(LINEAGE).stdout | head -n7 | awk 'BEGIN{printf "%s", "$*";}NR>1{printf "\t%s", $$2;}END{print "";}' > $@; \
	cp run_$*_BUSCOv3_proteome_$(LINEAGE)/short_summary_$*_BUSCOv3_proteome_$(LINEAGE).txt $(word 3,$^)/short_summary_$*.txt

%_proteome.fasta : %_CDSes.fasta logs
	@echo "Translating CDSes into proteins for $*"
	transeq -trim -sequence $(word 1,$^) -outseq $@ 2> $(word 2,$^)/transeq_$*.stderr

$(UNPUBCDSES) : %_unpub_CDSes.fasta : ../refs/Dsim_w501_v3_unpublished/%.fasta ../annotations/Dsim_w501_v3_unpublished/%.gff3 logs
	@echo "Extracting CDSes from $*_unpub"
	../tools/constructCDSesFromGFF3.pl -i $(word 1,$^) -g $(word 2,$^) -l 2> $(word 3,$^)/cCFG_$*_unpub.stderr > $@

$(EXTRACDSES) : %_CDSes.fasta : ../refs/old/%.fasta ../annotations/old/%.gff3 logs
	@echo "Extracting CDSes from $*"
	../tools/constructCDSesFromGFF3.pl -i $(word 1,$^) -g $(word 2,$^) -l 2> $(word 3,$^)/cCFG_$*.stderr > $@

$(CDSES) : %_CDSes.fasta : ../refs/%.fasta ../annotations/%.gff3 logs
	@echo "Extracting CDSes from $*"
	../tools/constructCDSesFromGFF3.pl -i $(word 1,$^) -g $(word 2,$^) -l 2> $(word 3,$^)/cCFG_$*.stderr > $@

$(SUBDIRS) :
	mkdir -p $@

clean :
	rm -f $(BUSCOS) $(EXTRABUSCOS) $(UNPUBBUSCOS)
	rm -f $(TRANBUSCOS) $(EXTRATRANBUSCOS) $(UNPUBTRANBUSCOS)
	rm -f $(PROTBUSCOS) $(EXTRAPROTBUSCOS) $(UNPUBPROTBUSCOS)
	rm -f $(CDSES) $(EXTRACDSES) $(UNPUBCDSES) *_proteome.fasta
	rm -f genome transcriptome proteome
	rm -rf BUSCO*summaries/ run_*/ tmp/ tmp_opt_*/ logs/
