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
#BUSCO result files:
BUSCOS := $(addsuffix _buscos.tsv,$(REFS))
#BUSCO result file(s) for old references (not our PacBio ones):
EXTRABUSCOS := $(addsuffix _buscos.tsv,$(OLDREFS))
#CDS FASTAs:
CDSES := $(addsuffix _CDSes.fasta,$(REFS))
#CDS FASTA(s) for old references (not our PacBio ones):
EXTRACDSES := $(addsuffix _CDSes.fasta,$(OLDREFS))
#Transcriptome BUSCO result files:
TRANBUSCOS := $(addsuffix _txome_buscos.tsv,$(REFS))
#Transcriptome BUSCO result file(s) for old references (not our PacBio ones):
EXTRATRANBUSCOS := $(addsuffix _txome_buscos.tsv,$(OLDREFS))
#Proteome BUSCO result files:
PROTBUSCOS := $(addsuffix _proteome_buscos.tsv,$(REFS))
#Proteome BUSCO result file(s) for old references (not our PacBio ones):
EXTRAPROTBUSCOS := $(addsuffix _proteome_buscos.tsv,$(OLDREFS))

.PHONY : genome_plot transcriptome_plot proteome_plot clean usage

.SECONDARY : $(BUSCOS) $(EXTRABUSCOS) genome transcriptome proteome

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
	@echo "clean -> Clean up all output files and directories"

genome_plot : genome
	generate_plot.py -wd BUSCO_summaries

genome : $(BUSCOS) $(EXTRABUSCOS)
	cat $^ > $@

#Run BUSCO v3 on each assembly, using 8 cores per job:
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

$(EXTRACDSES) : %_CDSes.fasta : ../refs/old/%.fasta ../annotations/old/%.gff3 logs
	@echo "Extracting CDSes from $*"
	constructCDSesFromGFF3.pl -i $(word 1,$^) -g $(word 2,$^) -l 2> $(word 3,$^)/cCFG_$*.stderr > $@

$(CDSES) : %_CDSes.fasta : ../refs/%.fasta ../annotations/%.gff3 logs
	@echo "Extracting CDSes from $*"
	constructCDSesFromGFF3.pl -i $(word 1,$^) -g $(word 2,$^) -l 2> $(word 3,$^)/cCFG_$*.stderr > $@

$(SUBDIRS) :
	mkdir -p $@

clean :
	rm -f $(BUSCOS) $(EXTRABUSCOS)
	rm -f $(TRANBUSCOS) $(EXTRATRANBUSCOS)
	rm -f $(PROTBUSCOS) $(EXTRAPROTBUSCOS)
	rm -f $(CDSES) $(EXTRACDSES)
	rm -f genome transcriptome proteome
	rm -rf BUSCO*summaries/ run_*/ tmp/ tmp_opt_*/ logs/
