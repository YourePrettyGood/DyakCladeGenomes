# DyakCladeGenomes

This repository contains all the scripts (in the `tools` subdirectory) and data (when size permits) used for analyzing constraint in the *Drosophila yakuba* clade based on PacBio genomes of *D. santomea*, *D. teissieri*, and *D. yakuba*.

## Dependencies common aross all segments:

1. [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit) is used implicitly in a lot of these pipelines, in particular the `fasta_formatter` tool.
1. Awk (GNU awk is preferred, haven't tested with POSIX awk)
1. Perl (used 5.26.2, should generally be compatible with 5.12+)
1. GNU CoreUtils (primarily `sort`, but also some others)
1. GNU Make (used 4.2.1)

## PacBio genomes:

We generated 4 reference-quality genome assemblies:

1. *D. santomea* STO CAGO 1482 (PacBio data from after 14 generations of inbreeding)
1. *D. teissieri* GT53w (PacBio data from after 14 generations of inbreeding)
1. *D. yakuba* NY73PB (NY73PB most closely resembles NY73 from Rogers et al. (2014) MBE, but has uncertain provenance)
1. *D. yakuba* Tai18E2 (Drosophila Species Stock Center #14021.0261-01)

Each line had over 70x of subreads, mostly generated from P6-C4 chemistry on a PacBio RS II at the UC-Irvine GHTF core (*D. teissieri* has some P5-C3 data).

Draft assemblies were generated using FALCON v0.7.0 and Canu 1.3, then these assemblies were merged with [quickmerge](https://github.com/mahulchak/quickmerge/) using the FALCON assembly as reference, and further scaffolding, gap-filling, and repeat resolution was performed using [FinisherSC](https://github.com/kakitone/finishingTool/) (Using the Imoteph fork). Assemblies were polished with Quiver, and then with Pilon. As the *D. teissieri* assembly was the most fragmented, in collaboration with Russ Corbett-Detig's group, we generated Hi-C data for *D. teissieri* GT53w (~100x) and *D. yakuba* NY73PB (~10x).

Initially, we mapped these Hi-C reads using the [Arima Genomics Hi-C mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline/), and scaffolded the genome using [SALSA2](https://github.com/machinegun/SALSA/). However, several lines of evidence suggested SALSA produced several misjoins, the most major being joins across all chromosome arms out of order and orientation. Thus, we performed misassembly detection and scaffolding into chromosome arms as follows:

Polished contig assemblies were aligned to the *D. melanogaster* ISO1 reference (FlyBase release 6.26), as well as the *D. yakuba* Tai18E2 reference (FlyBase release 1.05) to further determine joins necessary for chromosome-arm-length scaffolds by synteny, taking into account known inversions between species.

Further misassembly detection and correction was performed based on cytological banding patterns projected from *D. melanogaster* to contigs based on single-copy orthogroups assigned based on draft annotations. These cytological banding patterns were confirmed for *D. teissieri* GT53w and *D. yakuba* NY73PB based on Hi-C data via contact maps using HiCExplorer.

Contigs were joined using `manualScaffold.pl` with 100 bp gaps (per NCBI requirements) based on AGPs found under the `contigs` subdirectory of `refs`.

## Annotations:

We generated CDS annotations for the final chromosome-arm-scale genomes using [BRAKER2](https://github.com/Gaius-Augustus/BRAKER/) (version 2.1.4, commit 50e76af), using species-specific RNAseq data mapped using STAR 2.5.4b in two-pass mode, as well as the full proteome of *D. melanogaster* (FlyBase release 6.26, excluding proteins containing U or X).

Work is ongoing to refine transcript models using [StringTie](https://github.com/gpertea/stringtie/), and potentially more models and UTRs based on transcripts assembled using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/). We may also be testing out the [Comparative Annotation Toolkit](https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit) in Augustus-TM(R) and Augustus-CGP modes.

## Current data freeze availability:

The current assembly and annotation freeze (2019/12/11) is available (please be courteous about our publication intent) on request from [Peter Andolfatto](https://andolfattolab.com).

Pseudoreferences of the population samples may be linked from here in the future (please e-mail for more details).

## Assembly contiguity assessments:

Dependencies:

1. [R](https://www.r-project.org) (used 3.5.1)
1. [tidyverse](https://www.tidyverse.org) (R package, used 1.2.1)

Assembly contiguity statistics are produced in tabular (TSV) form in the `ctg_lengths` subdirectory using the command:
```
make -f asm_summaries.mk all
```

The results table is called `DyakClade_asm_summary_stats.tsv`, although this name can be adjusted by editing the `SUMMARYFILE` variable near the top of the `asm_summaries.mk` file.

The contig length plot input can be generated using the command:
```
make -f ctg_lengths.mk all
```

## BUSCO completeness assessments:

Dependencies:

1. [BUSCO](https://busco.ezlab.org) (v3, I used v3.0.2)
1. [EMBOSS](http://emboss.sourceforge.net) (I used 6.6.0-8 installed via dnf)

The assemblies and annotations (coding transcriptome and proteome) are assessed for completeness using `assembly_BUSCO.mk` in the BUSCO subdirectory. This Makefile depends on BUSCO v3.0.2 being installed and the `run_BUSCO.py` and `generate_plot.py` scripts being accessible via your PATH. This is most easily arranged using a module. Generating the proteomes for BUSCO assessment also requires EMBOSS transeq.

All assemblies are assessed against the narrowest BUSCO lineage containing Drosophila: Diptera. For the manuscript, we used NCBI-BLAST 2.2.31+, HMMER 3.1b2, and Augustus 3.3.2 (commit 8b1b14a7489e4545e89c8725dc33268f6c2a9117), with R 3.5.1 and EMBOSS transeq 6.6.0.0.

Command used:

```
module load BUSCO/3.0.2
make -j5 -f assembly_BUSCO.mk genome_plot
make -j5 -f assembly_BUSCO.mk transcriptome_plot
#make -j5 -f assembly_BUSCO.mk proteome_plot
```

Note: The plot for genome_plot had to be remade (so generate_plot.py was rerun) due to slight modifications made to the output R script:
Within the `annotate()` call, I adjusted `y=3` to `y=99` and `hjust=1` to `hjust=0`. This made the BUSCO summary text right-justified within each bar, instead of the default left-justified.

## Population Genetic Data:

Population samples used here include:

*D. santomea*:
1. 34 "G1" females (first-generation progeny of wild-caught individuals) (Illumina PE150)
1. 16 "Synthetic Diploid" females (progeny of an outbreeding cross of two isofemale lines) (Illumina PE215)

*D. teissieri*:
1. 8 "Synthetic Diploid" females (Illumina PE215)
1. 5 further "Synthetic Diploid" females (Illumina PE150)

*D. yakuba*:
1. 8 "Synthetic Diploid" females from lines collected on Sao Tome (Illumina PE215)
1. 10 "Synthetic Diploid" females from lines collected in mainland Africa (Illumina PE150)
1. 21 females from isofemale lines (treated as haploid -- one allele is chosen at random for any heterozygous sites) (Illumina PE50, PE75, PE100, from Rogers et al. (2014) MBE)
1. 12 further "Synthetic Diploid" females from lines collected in mainland Africa (Illumina PE150) (some overlapping with the 21 isofemale lines)

After filtering for inversion state in Dmel and line overlap in Dyak, we ended up with haploid sample sizes of 82 for *D. melanogaster*, 20 for *D. simulans*, 100 for *D. santomea*, 26 for *D. teissieri*, and 66 for *D. yakuba*, totalling 295 (excluding the 6 references).

Variant calls were generated and filtered using the [Pseudoreference Pipeline](https://github.com/YourePrettyGood/PseudoreferencePipeline/), in particular mapping with [BWA-MEM](https://github.com/lh3/bwa/) 0.7.17-r1188, marking duplicates with Picard MarkDuplicates commit 1ccb775, realigning around indels with GATK IndelRealigner 3.4-46-gbc02625, calling variants with [BCFtools](https://github.com/samtools/bcftools/) 1.7-1-ge07034a (htslib 1.7-23-g9f9631d), and masking sites for which the following compound criterion held true:

`DP <= [half of post-markdup depth] || MQ <= 20.0 || QUAL <= 26.0`

Indels were excluded from analysis and generation of pseudoreferences.

## Determination of Single-copy Orthogroups:

Dependencies:

1. [OrthoFinder](https://github.com/davidemms/OrthoFinder) 2.3.8 (installed via bioconda)
1. Diamond 0.9.29 (installed as a dependency of OrthoFinder 2.3.8)
1. IQ-Tree 1.5.5 (installed as a dependency of OrthoFinder 2.3.8)
1. [Clustal Omega](http://clustal.org/omega) 1.2.4

Using the PacBio genomes and annotations, we established orthogroups (and an estimate of the species tree) using [OrthoFinder](https://github.com/davidemms/OrthoFinder/), using [Diamond](https://github.com/bbuchfink/diamond/) for all-vs-all comparisons, constructing orthogroup trees by aligning with [Clustal Omega](http://www.clustal.org/omega/), and inferring trees with [IQ-Tree](https://www.iqtree.org/). The `config.json` of OrthoFinder needed to be modified to support Clustal Omega.

Only the longest CDS for each gene was extracted and translated for input into OrthoFinder for each line. Length ties were broken by lexicographical order of the transcript ID to ensure consistency between runs.

Single-copy orthogroups were established simply by extracting only those orthogroups with a single member transcript from each input proteome.

All of this is achieved by running the following commands while in the `OrthoFinder_clustalo` subdirectory:

```
#Using /usr/bin/time -v to get runtime and memory stats, and redirecting
# STDERR and STDOUT to avoid terminal spam
#Prep everything for OrthoFinder:
#Took 
/usr/bin/time -v make -j4 -f prep_proteomes.mk proteome 2> prep_proteomes_proteome_j4.stderr > prep_proteomes_proteome_j4.stdout

#Run OrthoFinder:
#Installed into its own conda environment
#For reference, it took about 6 minutes to run Diamond and generate Orthogroups.tsv
#However, it took about 17 hours to completely finish (i.e. all gene
# tree inference and the final species tree inference with STRIDE)
#In principle, you only need to wait for Orthogroups.tsv to be made to
# continue with the rest of this analysis.
#Took 16h 42m 26s, 583% CPU, 8,116,956 KB Max RSS
#Spent 12h 45m 29s on ModelFinder testing 546 models on the concatenated
# SpeciesTree alignment -- this is the rate-limiting step
source activate OrthoFinderenv
/usr/bin/time -v orthofinder -t 36 -a 1 -M msa -S diamond -A clustalo -T iqtree -f OrthoFinder_proteomes 2> logs/orthofinder_diamond_clustalo_IQTree_MSSTYY_t36.stderr > logs/orthofinder_diamond_clustalo_IQTree_MSSTYY_t36.stdout

#Remap the Orthogroups.tsv file:
/usr/bin/time -v make -f prep_proteomes.mk SCOmap 2> logs/prep_proteomes_SCOmap.stderr > logs/prep_proteomes_SCOmap.stdout

#Do a final check of SCO count per OrthoFinder, per awk, and the overlap of the two:
#The three numbers in the second column of the SCOcheck file should all match
/usr/bin/time -v make -f prep_proteomes.mk SCOcheck 2> logs/prep_proteomes_SCOcheck.stderr > logs/prep_proteomes_SCOcheck.stdout
```

## CDS alignments:

Dependencies:

1. [PRANK](https://wasabiapp.org/software/prank/) v.170427

The longest CDS for each gene per species was extracted, and these were parsed into FASTAs by orthogroup -- only single-copy orthogroups were preserved. These CDSes were then trimmed for stop codons, and aligned with [PRANK](https://wasabiapp.org/software/prank/). These same CDSes were extracted and parsed from the pseudoreferences, where any diploid samples had their sequences randomly split into two haplotypes, and were trimmed for stop codons. Gaps were inserted into these pseudoreference sequences based on the gaps inserted by PRANK into their corresponding species reference CDS, generating an alignment of all the references and their corresponding pseudoreferences. Since the pseudoreferences only updated SNPs, it was unnecessary (and computationally inefficient) to include them in the PRANK alignment. Due to the way the stop codon trimming worked, we had to further square-off the alignment after adding the pseudoreferences, as some pseudoreferences would have the stop codon masked (hence not being recognized by the trimming script).

All of this is achieved by running the following commands while in the `CDS` subdirectory:

```
#Using /usr/bin/time -v to get runtime and memory stats, and redirecting
# STDERR and STDOUT to avoid terminal spam
#-j16 is to parallelize across 16 threads (only 5 references to deal with)
#This parses out the reference CDSes:
#Took 0m 52.31s, 115% CPU, 163,244 KB Max RSS
/usr/bin/time -v make -j16 -f align_CDSes.mk parse_ref_CDSes 2> parse_ref_CDSes_j16.stderr > parse_ref_CDSes_j16.stdout

#-j36 is to parallelize across 36 threads (server had 40, there are 197 pseudorefs)
#This parses out the pseudoref CDSes:
#Took 5m 11.12s, 2270% CPU, 148,948 KB Max RSS
/usr/bin/time -v make -j36 -f align_CDSes.mk parse_sample_CDSes 2> parse_sample_CDSes_j36.stderr > parse_sample_CDSes_j36.stdout

#Again with -j36, but now we have >9k alignments to run
#This aligns all the 9,375 reference orthogroups (slowest part), and
# adds in the pseudoref haplotypes after (fast):
#Took 5h 51m 51s, 3129% CPU, 531,512 KB Max RSS
/usr/bin/time -v make -j36 -f align_CDSes.mk align_all 2> align_all_j36.stderr > align_all_j36.stdout
```

These alignments were further filtered downstream for missing data.

## Intron alignments:

Dependencies:

1. [Clustal Omega](http://clustal.org/omega) 1.2.4

All introns corresponding to transcripts in the single-copy orthogroup file were extracted from each reference. During extraction, the header of a given intron also contains information about the length of the intron, and the lengths of its flanking exons. We use a threshold of a maximum 20% relative difference between flanking exon lengths to establish that two introns are homologous (along with belonging to members of the same orthogroup). We then only keep intron groups with a single matching orthologous intron for each species. The reference sequences of these orthologous intron groups are parsed into intron-specific FASTAs, and aligned with Clustal Omega. As with the CDSes, after adjusting the intron map to only include species ID and not line ID, we extract and parse pseudoreference introns (again, splitting into pseudohaplotypes randomly), we add gaps according to species based on the reference intron alignments.

All of this is achieved by running the following commands while in the `introns` subdirectory:

```
#Using /usr/bin/time -v to get runtime and memory stats, and redirecting
# STDERR and STDOUT to avoid terminal spam
#-j16 is to parallelize across 16 threads (only 5 references to deal with)
#This parses out the reference CDSes:
#Took 2m 36.35s, 109% CPU, 148,848 KB Max RSS
/usr/bin/time -v make -j16 -f align_introns.mk parse_ref_introns 2> parse_ref_introns_j16.stderr > parse_ref_introns_j16.stdout

#-j36 is to parallelize across 36 threads (server had 40, there are 197 pseudorefs)
#This parses out the pseudoref introns:
#Took 8m 56.47s, 1260% CPU, 148,252 KB Max RSS
/usr/bin/time -v make -j36 -f align_introns.mk parse_sample_introns 2> parse_sample_introns_j36.stderr > parse_sample_introns_j36.stdout

#Again with -j36, but now we have 18,444 alignments to run
#This aligns all the orthologous intron groups (slowest part), and
# adds in the pseudoref haplotypes after (fast):
#Took 2h 38m 24s, 3663% CPU, 56,230,532 KB Max RSS
/usr/bin/time -v make -j36 -f align_introns.mk align_all 2> align_all_j36.stderr > align_all_j36.stdout
```

These alignments were further filtered downstream for missing data.

## Synteny dotplots

Dependencies:

1. [MUMmer](https://github.com/mummer4/mummer) (You could use MUMmer 3.23 if you want, without any issues.)
1. GNUplot (MUMmer depends on this for mummerplot functionality)



## Genome size estimation

Dependencies:

1. `calculate_read_depth_gc_windows.py` and `adjust_read_depth_windows.py` from [John Davey's scripts from Davey et al. (2016) G3](https://github.com/johnomics/Heliconius_melpomene_version_2), with some modifications to catch divide-by-zero exceptions
1. [BEDtools](https://github.com/arq5x/bedtools2) (used v2.27.1-1-gb87c4653)
1. [Pseudoreference Pipeline](https://github.com/YourePrettyGood/PseudoreferencePipeline)
