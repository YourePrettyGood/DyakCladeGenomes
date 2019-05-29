# DyakCladeGenomes

This repository contains all the scripts (in the `tools` subdirectory) and data (when size permits) used for analyzing constraint in the *Drosophila yakuba* clade based on PacBio genomes of *D. santomea*, *D. teissieri*, and *D. yakuba*.

## PacBio genomes:

We generated 4 reference-quality genome assemblies:

1. *D. santomea* STO CAGO 1482 (PacBio data from after 14 generations of inbreeding)
1. *D. teissieri* GT53w (PacBio data from after 14 generations of inbreeding)
1. *D. yakuba* NY73PB (NY73PB most closely resembles NY73 from Rogers et al. (2014) MBE, but has uncertain provenance)
1. *D. yakuba* Tai18E2 (Drosophila Species Stock Center #14021.0261-01)

Each line had over 70x of subreads, mostly generated from P6-C4 chemistry on a PacBio RS II at the UC-Irvine GHTF core (*D. teissieri* has some P5-C3 data).

Draft assemblies were generated using FALCON v0.7.0 and Canu 1.3, then these assemblies were merged with [quickmerge](https://github.com/mahulchak/quickmerge/) using the FALCON assembly as reference, and further scaffolding, gap-filling, and repeat resolution was performed using [FinisherSC](https://github.com/kakitone/finishingTool/) (Using the Imoteph fork). Assemblies were polished with Quiver, and then with Pilon. As the *D. teissieri* assembly was the most fragmented, in collaboration with Russ Corbett-Detig's group, we generated Hi-C data for *D. teissieri* GT53w (~100x) and *D. yakuba* NY73PB (~10x). We mapped these Hi-C reads using the [Arima Genomics Hi-C mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline/), and scaffolded the genome using [SALSA2](https://github.com/machinegun/SALSA/).

Polished contig assemblies were aligned to the *D. melanogaster* ISO1 reference (FlyBase release 6.26), as well as the *D. yakuba* Tai18E2 reference (FlyBase release 1.06) to further determine joins necessary for chromosome-arm-length scaffolds by synteny, taking into account known inversions between species. Contigs were joined using `manualScaffold.pl`.

## Annotations:

We generated CDS annotations for the final chromosome-arm-scale genomes using [BRAKER2](https://github.com/Gaius-Augustus/BRAKER/) (version 2.1.0, with slight modifications), using species-specific RNAseq data mapped using STAR 2.5.2a in two-pass mode, as well as the full proteome of *D. melanogaster* (FlyBase release 6.26).

Work is ongoing to refine transcript models using [StringTie](https://github.com/gpertea/stringtie/), and potentially more models and UTRs based on transcripts assembled using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/).

## Population Genetic Data:

Population samples used here include:

*D. santomea*:
1. 34 "G1" females (first-generation progeny of wild-caught individuals) (Illumina PE150)
1. 16 "Synthetic Diploid" females (progeny of an outbreeding cross of two isofemale lines) (Illumina PE215)

*D. teissieri*:
1. 8 "Synthetic Diploid" females (Illumina PE215)

*D. yakuba*:
1. 8 "Synthetic Diploid" females from lines collected on Sao Tome (Illumina PE215)
1. 10 "Synthetic Diploid" females from lines collected in mainland Africa (Illumina PE150)
1. 21 females from isofemale lines (treated as haploid -- one allele is chosen at random for any heterozygous sites) (Illumina PE50, PE75, PE100, from Rogers et al. (2014) MBE)

Variant calls were generated and filtered using the [Pseudoreference Pipeline](https://github.com/YourePrettyGood/PseudoreferencePipeline/), in particular mapping with [BWA-MEM](https://github.com/lh3/bwa/), marking duplicates with Picard MarkDuplicates, realigning around indels with GATK IndelRealigner, calling variants with [BCFtools](https://github.com/samtools/bcftools/), and masking sites for which the following compound criterion held true:

`DP <= [half of post-markdup depth] || MQ <= 20.0 || QUAL <= 26.0`

Indels were excluded from analysis and generation of pseudoreferences.

## Determination of Single-copy Orthogroups:

Using the PacBio genomes and annotations, we established orthogroups (and an estimate of the species tree) using [OrthoFinder](https://github.com/davidemms/OrthoFinder/) v2.2.7 (via Bioconda), using [Diamond](https://github.com/bbuchfink/diamond/) 0.9.21 for all-vs-all comparisons, constructing orthogroup trees by aligning with [Clustal Omega](https://www.clustal.org/omega/) 1.2.4, and inferring trees with [IQ-Tree](https://www.iqtree.org/) 1.5.5. The `config.json` of OrthoFinder needed to be modified to support Clustal Omega.

Only the longest CDS for each gene was extracted and translated for input into OrthoFinder for each line. Length ties were broken by lexicographical order of the transcript ID to ensure consistency between runs.

Single-copy orthogroups were established simply by extracting only those orthogroups with a single member transcript from each input proteome.

All of this is achieved by running the following commands while in the `OrthoFinder_clustalo` subdirectory:

```
#Prep everything for OrthoFinder:
make -f prep_proteomes.mk proteome

#Run OrthoFinder:
#Installed into its own conda environment
#For reference, it took about 1-2 hours to run Diamond and generate Orthogroups.csv
#However, it took a bit under 12 hours to completely finish (i.e. all gene
# tree inference and the final species tree inference with STRIDE)
#In principle, you only need to wait for Orthogroups.csv to be made to
# continue with the rest of this analysis.
#Took 11h 38m 39s, 643% CPU, 7,256,424 KB Max RSS
source activate OrthoFinderenv
nohup /usr/bin/time -v orthofinder -t 32 -a 1 -M msa -S diamond -A clustalo -T iqtree -f OrthoFinder_proteomes 2> logs/orthofinder_diamond_clustalo_IQTree_MSTYY.stderr > logs/orthofinder_diamond_clustalo_IQTree_MSTYY.stdout

#Remap the Orthogroups.csv file:
make -f prep_proteomes.mk SCOmap

#Do a final check of SCO count per OrthoFinder, per awk, and the overlap of the two:
#The three numbers in the second column of the SCOcheck file should all match
make -f prep_proteomes.mk SCOcheck
```

## CDS alignments:

The longest CDS for each gene per species was extracted, and these were parsed into FASTAs by orthogroup -- only single-copy orthogroups were preserved. These CDSes were then trimmed for stop codons, and aligned with [PRANK](https://wasabiapp.org/software/prank/) v.170427. These same CDSes were extracted and parsed from the pseudoreferences, where any diploid samples had their sequences randomly split into two haplotypes, and were trimmed for stop codons. Gaps were inserted into these pseudoreference sequences based on the gaps inserted by PRANK into their corresponding species reference CDS, generating an alignment of all the references and their corresponding pseudoreferences. Since the pseudoreferences only updated SNPs, it was unnecessary (and computationally inefficient) to include them in the PRANK alignment. Due to the way the stop codon trimming worked, we had to further square-off the alignment after adding the pseudoreferences, as some pseudoreferences would have the stop codon masked (hence not being recognized by the trimming script).

All of this is achieved by running the following commands while in the `CDS` subdirectory:

```
#Using /usr/bin/time -v to get runtime and memory stats, and redirecting
# STDERR and STDOUT to avoid terminal spam
#-j4 is to parallelize across 4 threads (only 4 references to deal with)
#This parses out the reference CDSes:
#Took 1m 9.53s, 107% CPU, 148,032 KB Max RSS
/usr/bin/time -v make -j4 -f align_CDSes.mk parse_ref_CDSes 2> parse_ref_CDSes.stderr > parse_ref_CDSes.stdout

#-j30 is to parallelize across 30 threads (server had 40, there are 97 pseudorefs)
#This parses out the pseudoref CDSes:
#Took 46.59s, 765% CPU, 57,444 KB Max RSS
/usr/bin/time -v make -j30 -f align_CDSes.mk parse_sample_CDSes 2> parse_sample_CDSes.stderr > parse_sample_CDSes.stdout

#Again with -j30, but now we have >9k alignments to run
#This aligns all the 9,580 reference orthogroups (slowest part), and
# adds in the pseudoref haplotypes after (fast):
#Took 6h 24m 32s, 1688% CPU, 728,936 KB Max RSS
/usr/bin/time -v make -j30 -f align_CDSes.mk align_all 2> align_all.stderr > align_all.stdout
```

These alignments were further filtered downstream for missing data.

## Intron alignments:

All introns corresponding to transcripts in the single-copy orthogroup file were extracted from each reference. During extraction, the header of a given intron also contains information about the length of the intron, and the lengths of its flanking exons. We use a threshold of a maximum 20% relative difference between flanking exon lengths to establish that two introns are homologous (along with belonging to members of the same orthogroup). We then only keep intron groups with a single matching orthologous intron for each species. The reference sequences of these orthologous intron groups are parsed into intron-specific FASTAs, and aligned with Clustal Omega. As with the CDSes, after adjusting the intron map to only include species ID and not line ID, we extract and parse pseudoreference introns (again, splitting into pseudohaplotypes randomly), we add gaps according to species based on the reference intron alignments.

All of this is achieved by running the following commands while in the `introns` subdirectory:

```
#Using /usr/bin/time -v to get runtime and memory stats, and redirecting
# STDERR and STDOUT to avoid terminal spam
#-j10 is to parallelize across 10 threads (only 4 references to deal with)
#This parses out the reference CDSes:
#Took 2m 40.47s, 105% CPU, 142,704 KB Max RSS
/usr/bin/time -v make -j10 -f align_introns.mk parse_ref_introns 2> parse_ref_introns.stderr > parse_ref_introns.stdout

#-j10 is to parallelize across 10 threads (server had 40, there are 97 pseudorefs)
#This parses out the pseudoref introns:
#Took 2m 1.41s, 285% CPU, 56,240 KB Max RSS
/usr/bin/time -v make -j10 -f align_introns.mk parse_sample_introns 2> parse_sample_introns.stderr > parse_sample_introns.stdout

#Again with -j32, but now we have 19,978 alignments to run
#This aligns all the orthologous intron groups (slowest part), and
# adds in the pseudoref haplotypes after (fast):
#Took 2h 0m 50s, 3611% CPU, 36,500,264 KB Max RSS
/usr/bin/time -v make -j32 -f align_introns.mk align_all 2> align_all.stderr > align_all.stdout
```

These alignments were further filtered downstream for missing data.
