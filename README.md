# DyakCladeGenomes

These are the scripts and programs that were used to generate, annotate, process, and analyze PacBio genomes (and popgen samples mapped to them) of the species of the *Drosophila yakuba* clade, i.e. *D. santomea*, *D. teissieri*, and *D. yakuba* (NY73PB line was used as reference for Dyak popgen samples).

The Perl and C++ scripts that I've written have a `-h` option, so take a look, and let me know if anything is unclear! (Awk scripts here don't have a help flag.)

## Categories
1. [De novo assembly-related](README.md#de-novo-assembly-related-scripts)
1. [Annotation-related](README.md#annotation-related-scripts)
1. [Genome-wide statistics](README.md#genome-wide-statistics-programsscripts)
1. [MSA-related](README.md#msa-related-scripts)

## Dependencies in example usages:

All of these scripts and programs are standalone, but in typical usage they go hand in hand with some very useful external tools:

1. [FASTX-Toolkit](https://hannonlab.cshl.edu/fastx_toolkit/) from the Hannon Lab at CSHL
1. [GenomeTools](https://github.com/genometools/genometools/)
1. [BEDtools](https://github.com/arq5x/bedtools2/)
1. [SAMtools](https://github.com/samtools/samtools/)
1. Many of these tools are designed to work on pseudoreferences produced by the [Pseudoreference Pipeline](https://github.com/YourePrettyGood/PseudoreferencePipeline/)
1. [GNU Parallel](https://www.gnu.org/software/parallel/) is used to parallelize many of these tasks that operate on the thousands of orthologs and tens of thousands of introns, etc.

## One-liners:

### Convert FASTQ to unwrapped FASTA (FASTQ can be gzipped, or uncompressed):

`gzip -dcf [FASTQ file] | awk 'NR%4==1||NR%4==2' | sed 's/@/>/' > [FASTA file]`

### Output lengths of longest `N` contigs or scaffolds from an assembly:

`[path to FASTX toolkit]/fasta_formatter -i [assembly FASTA] | awk '/^>/{scafname=$1;}!/^>/{print substr(scafname, 2)"\t"length($0);}' | sort -k2,2nr | head -n[N]`

**Note: FASTX toolkit is available from [Hannon Lab at CSHL](http://hannonlab.cshl.edu/fastx_toolkit/)**

### Split scaffolds into contig at N gaps of length >= 1:

`[path to FASTX toolkit]/fasta_formatter -i [assembly FASTA] | awk '/^>/{header=$0;}!/^>/{split($0, ctgs, /[N]+/); for (ctg in ctgs) {print header"_ctg"ctg; print ctgs[ctg];};}' > [unwrapped contigs FASTA]`

Just adjust the regex in the `split()` call if you want to split at gaps of specified minimum length.

### Generate a BED of gaps from a genome FASTA:

`[path to FASTX toolkit]/fasta_formatter -i [assembly FASTA] | perl -e 'use warnings; use strict; my $scafname = ""; while (my $line = <>) { chomp $line; if ($line =~ /^>/) {$scafname = substr $line, 1;} else { my $seq = $line; print join("\t", $scafname, $-[0], $+[0]+1), "\n" while $seq =~ /[^ACGTacgt]+/g;};}' > [gaps BED file]`

This outputs BED intervals for any regions with characters other than ACGT, so degenerate bases will be included in these intervals. Just adjust the regex if you want different gap detection behaviour (e.g. `/[Nn]+/g` if you only want to consider N gaps). The regex I used above was intended for use with a haploid reference assembly for the genome size estimation pipeline from Davey et al. (2016) G3 [doi:10.1534/g3.115.023655](https://dx.doi.org/10.1534/g3.115.023655)

## De novo assembly-related scripts:

### `NX.pl`

Calculates some standard contiguity statistics from an assembly FASTA.  By default, it calculates the N50 and L50 (length-weighted median contig/scaffold size and number, respectively), and all calls to it also indicate the total assembly size, number of contigs/scaffolds in the assembly, length of longest and shortest contig/scaffold, and average contig/scaffold length.  Options include changing the quantile from 50 (e.g. set quantile to 90 to calculate N90 and L90, setting it to 50 is the same as default), and specifying the genome size (e.g. in the default case, to calculate an NG50).

Multiple quantiles can be evaluated by specifying a comma-separated (no spaces) list of quantiles as input.

You can pass `-` as the FASTA file path, and the script will read the FASTA from STDIN, so it can easily be incorporated into piped one-liners.

Also, interestingly, with some upstream fiddling to get reads into FASTA format, you can calculate the NRX statistics (like NR25 and NR40) with this script, although it's fairly slow for large numbers of reads.  We're limited by Perl's sort routine's time complexity, unfortunately.

Example usage:

Calculate N50, L50, and other stats for Drosophila melanogaster FlyBase release 6.13 assembly:

`NX.pl dmel-all-chromosome-r6.13.fasta`

Calculate N90, L90, etc. for the same assembly:

`NX.pl dmel-all-chromosome-r6.13.fasta 90`

Calculate the NG50, LG50, etc., assuming G=130,000,000 bp:

`NX.pl dmel-all-chromosome-r6.13.fasta 50 130000000`

Calculate N50 on a gzipped assembly (e.g. from 10x Genomics' Supernova mkoutput):

`gzip -dc MySupernovaAsm.fa.gz | NX.pl -`

or

`NX.pl <(gzip -dc MySupernovaAsm.fa.gz)`

or

`NX.pl MySupernovaAsm.fa.gz`

though piping and command substitution seem to be faster in most cases.

Calculate the N10, N50, and N90 for the Dmel assembly used above:

`NX.pl dmel-all-chromosome-r6.13.fasta 10,50,90`

### `manualScaffold.pl`

Usage:

`manualScaffold.pl -i [path to unscaffolded FASTA] [options] <configuration string>`

This was originally developed as a quick way to manually scaffold contigs into chromosome arms.  You must know *a priori* what the order and orientation of contigs needs to be, as you specify a configuration string, or supply an AGP file to dictate how the script sews the contigs together, and how to label the resultant scaffolds.

If you wish to read the input unscaffolded FASTA from STDIN, just omit the `-i` option.

The `-u` option prints contigs that were not scaffolded into the output FASTA as well.  It just saves you from having to write out a very long configuration string or AGP file.

It has now been adapted to scaffold based on an [AGP version 2.0](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/) file with the `-a` option, which properly treats `N` and `U` gap records for length, although this script ignores columns 2, 3, 7, and 8 (so it takes the entirety of the contig even if the AGP says to take a substring).

The original configuration string format isn't the most intuitive, but it was easy to come up with.  Defined using the [augmented BNF grammar](https://en.wikipedia.org/wiki/Augmented_Backus%E2%80%93Naur_form):

`config-string = scaffold*(scaffold-delim scaffold)`

`scaffold-delim = "=>"`

`scaffold = scaffold-name record-delim contig*(contig-delim contig)`

`record-delim = ":"`

`contig-delim = "->"`

`contig = contig-name [revcomp]`

`revcomp = "*"`

`contig-name = *(ALPHA / DIGIT / "_" / "|")`

So for example, if I had 3 contigs generated by Quiver polishing of FinisherSC-improved contigs, and I knew that they should be given the scaffold name chr1 and joined as Segkk0|quiver's end with Segkk1|quiver's beginning, and Segkk1|quiver's end with Segkk2|quiver's end, the configuration string would be:

`chr1:Segkk0|quiver->Segkk1|quiver->Segkk2|quiver*`

If there were more scaffolds to make, the configuration string might look like:

`chr1:Segkk0|quiver->Segkk1|quiver->Segkk2|quiver*=>chr2:Segkk3|quiver*->Segkk19|quiver->Segkk8|quiver*`

In fact, the configuration string I used for *D. yakuba* Tai18E2 before Pilon polishing was:

`X:Segkk0_quiver*->Segkk61_quiver*->Segkk60_quiver->Segkk41_quiver*->Segkk49_quiver=>2L:Segkk46_quiver*->Segkk47_quiver=>2R:Segkk7_quiver*=>3L:Segkk1_quiver*->Segkk54_quiver=>3R:Segkk48_quiver->Segkk55_quiver=>4:Segkk3_quiver*`

### `configStringToAGP.pl`

Usage:

`configStringToAGP.pl -i [.fai file] [-c config string file] [-u] [-a output AGP] <config string>`

`-i` specifies the FASTA index (.fai) for the unscaffolded FASTA, which is necessary to detect contigs missing from the config string.

`-c` acts as an alternate to providing the config string as a positional argument.

`-u` ensures that contigs missing from the config string are included in the output AGP.

`-a` specifies the path for the output AGP. If omitted, this defaults to STDOUT.

As you can imagine, this script basically does one part of what `prepContigsForGenBank.pl` does.  In fact it was derived from the same code, I just didn't need a FASTQ, just an easy way to get an AGP for quick manual scaffolding.

Note the slight differences between an assembly scaffolded with an AGP produced by `configStringToAGP.pl` versus scaffolding with the config string: The AGP version 2.0 specification states that gaps of unknown size (`U` records) should have length 100, whereas `manualScaffold.pl` using a config string will by default make gaps of 500 Ns.

### `softmaskFromHardmask.cpp`

This program will softmask a genome, given an unmasked genome and a hardmasked genome. The underlying code is quite simple, but it turns out to be a lot easier to use this than to re-run RepeatMasker on softmasking mode if you accidentally ran hardmasking mode. Just make sure both the unmasked and hardmasked FASTAs are either fully unwrapped, or wrapped to exactly the same length, otherwise the program will error out. Also, it doesn't seem to handle process substitutions correctly (e.g. will error out about different line wrappings if you feed it process substitutions calling fasta_formatter).

Usage:

`softmaskFromHardmask [unmasked FASTA] [hardmasked FASTA] > [softmasked FASTA]`

## Annotation-related scripts:

### `CDStoGenomicIntervals.pl`

Usage:

`CDStoGenomicIntervals.pl -i [CDS-space BED] -g [GFF3]`

`-i` specifies a BED file of intervals in CDS-space which you want to convert to genome-space.

`-g` specifies the GFF3 file that provides the mapping from CDS-space to genome-space.

Example Usage:

`CDStoGenomicIntervals.pl -i Dyak_NY73PB_v2_4fold_sites.bed -g Dyak_NY73PB_v2_BRAKER2_RNAseq_DmelProteins_adjusted_sorted.gff | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > Dyak_NY73PB_v2_4fold_sites_genomic.bed`

This script takes a BED of CDS-space intervals and a GFF3 annotation (requiring CDS records that have a Parent= tag whose value would match the scaffold ID of the CDS-space BED), and outputs the equivalent BED in genome-space, accounting for intron-spanning intervals.

The output BED is **not** sorted or merged, so pass the output through `sort -k1,1 -k2,2n -k3,3n` and `bedtools merge -i -` before downstream usage.

### `codingSitesByDegeneracy.pl`

Usage:

`codingSitesByDegeneracy.pl -f [Degree of Degeneracy] -i [CDS FASTA]`

`-f` specifies the degree of degeneracy of sites you want to identify. This degree is the number of nucleotides for which, when a nucleotide is substituted at this position, the amino acid encoded by the codon remains the same as the input codon. For example, at a 4-fold site (degree of degeneracy 4), any of the 4 nucleotides may be substituted at this position to obtain the same amino acid.

`-i` specifies the path to a FASTA of coding sequences (CDSes) whose sites you want to categorize.

Example Usage:

`codingSitesByDegeneracy.pl -f 4 -i Dyak_NY73PB_v2_BRAKER2_RNAseq_DmelProteins_CDSes.fasta | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > Dyak_NY73PB_v2_4fold_sites.bed`

This script takes a degree of degeneracy (specified by the `-f` flag, e.g. `-f 4` for 4-fold synonymous sites, `-f 2` for 2-fold sites, etc.) and a FASTA of coding sequences (CDSes, say generated by `constructCDSesFromGFF3.pl`) and outputs an unsorted and unmerged BED indicating the CDS-space positions of sites with the specified degree of degeneracy.

When coupled with `CDStoGenomicIntervals.pl`, `subsetVCFstats.pl`, and `calculateDxy`, it allows for calculation of an uncorrected per-site dS (or even dN for `-f 0` or `-f 1`).

### `constructCDSesFromGFF3.pl`

Usage:

`constructCDSesFromGFF3.pl -i [FASTA] -g [GFF3] [-e]`

`-i` specifies the FASTA from which to extract sequences. `-i` is optional, defaults to STDIN

`-g` specifies the GFF3, which must be in the same coordinate space as the FASTA, and must be locally sorted (i.e. exons and CDSes within a gene must be in ascending positional order).

`-e` is optional, outputs an extra 2 lines per FASTA record, which are the "Exon Range String" and the "CDS Range String".  These indicate the coordinates of the boundaries of each CDS record in genomic-space and CDS-space, respectively. They are always written with the initial exon first, and terminal exon last, so CDSes on the negative strand have Exon Range String elements listed in reverse order, and with larger number first in each range.

`-l` Only extracts the longest isoform for each gene, and skips any remaining isoforms for that gene. This is particularly useful for use as input to ortholog-finding pipelines like [OrthoFinder](https://github.com/davidemms/OrthoFinder/).

This script takes in a GFF3 file (e.g. genome annotation -- be very certain that it is GFF3 and not GFF2 or GTF), and a genome FASTA, and outputs all CDSes from the annotation as unwrapped FASTA.

I've used this to prepare transcriptomes (minus UTRs) for differential expression analysis, and it does not suffer from the off-by-one error you get from bedtools getfasta (plus you don't get all the extraneous individual exon and CDS FASTA records).

For instance, if you have the "GFF" output from the BRAKER1 pipeline, this is not actually GFF, but an oddly formatted GTF, so you can convert from GTF to GFF3 using the Augustus script `gtf2gff.pl` as follows:

```
gtf2gff.pl --printExon --printUTR --gff3 --out BRAKER_output_gtf2gff.gff3 < BRAKER_output.gff
constructCDSesFromGFF3.pl -i BRAKER_genome.fasta -g BRAKER_output_gtf2gff.gff3 > BRAKER_output_gtf2gff_transcripts.fasta
```

Oftentimes with BRAKER1- or BRAKER2-derived GFF3 files, you'll want to fix and fully sort them. My typical command set for getting a valid, sorted GFF3 ready for this script is:

```
awk 'BEGIN{FS="\t";OFS="\t";}$2=="AUGUSTUS"{if ($3 == "transcript") {split($9, genetxarr, "."); $9="gene_id \""genetxarr[1]"\"; transcript_id \""$9"\";";}; print $0;}' braker/[BRAKER run ID]/augustus.hints.gtf | ~/bin/augustus-3.3/augustus/scripts/gtf2gff.pl --printExon --gff3 --out=[annotation ID].gff3 2> [annotation ID]_gtf2gff.stderr
samtools faidx [softmasked genome FASTA]
cat <(awk 'BEGIN{print "##gff-version 3";}{print "##sequence-region "$1" 1 "$2;}' [softmasked genome FASTA].fai) <(fillMissingGenes.awk [annotation ID].gff3) > [annotation ID]_adjusted.gff3
gt gff3 -sort -tidy -retainids -checkids [annotation ID]_adjusted.gff3 2> [annotation ID]_adjusted_sorting.stderr > [annotation ID]_adjusted_sorted.gff3
```

### `extractIntronsFromGFF3.pl`

Usage:

`extractIntronsFromGFF3.pl -i [FASTA] -g [GFF3]`

`-i` specifies the FASTA from which to extract sequences. `-i` is optional, defaults to STDIN

`-g` specifies the GFF3, which must be in the same coordinate space as the FASTA, and must be locally sorted (i.e. exons and CDSes within a gene must be in ascending positional order).

`-b` specifies which type of GFF3 feature to use as boundaries of the inferred introns. (Default: exon, but for consistency between annotations, you may want to use CDS, as some annotations may be missing UTR exons)

This script takes in a GFF3 file (e.g. genome annotation -- be very certain that it is GFF3 and not GFF2 or GTF), and a genome FASTA, and outputs all CDSes from the annotation as unwrapped FASTA.

The FASTA headers are space-separated lists of information useful for downstream filtering for homology: Intron ID (GFF3-style), intron coordinates (FlyBase-style), intron length, left flanking exon length, right flanking exon length, left flanking exon ID, right flanking exon ID.

Beware that the extracted intronic sequence will contain the two splice junctions, so you'll probably want to trim those before aligning.

GFF3-style intron and exon IDs are of the form: [transcript ID].intron# or [transcript ID].exon#

FlyBase-style coordinates are of the form: `[scaffold ID]:[left bound]..[right bound]([strand])`

### `fillMissingGenes.awk`

This awk script takes a GFF3 file produced by gtf2gff.pl (from Augustus, tested with Augustus version 3.3), and fills in any missing `gene` and `mRNA` records, which happens surprisingly often. This script has primarily been used in the process of converting the `augustus.hints.gtf` file produced by the BRAKER2 (v2.1.0) pipeline into a specification-compliant GFF3.

Example Usage:

`cat <(awk 'BEGIN{print "##gff-version 3";}{print "##sequence-region "$1" 1 "$2;}' Dsim_w501_Pilon_chromosome_arms_mtDNA_softmasked_w60.fasta.fai) <(fillMissingGenes.awk Dsim_w501_BRAKER2_RNAseq_DmelProteins.gff3) > Dsim_w501_BRAKER2_RNAseq_DmelProteins_adjusted.gff3`

(Note that the first awk command in the example is simply adding in the appropriate sequence-region header lines based on the FASTA index file, as well as the GFF version header.)

The resultant GFF3 (with GFF version and sequence-region headers) should now be compatible with GenomeTools' sorting algorithm, like so:

`gt gff3 -sort -tidy -retainids -checkids Dsim_w501_BRAKER2_RNAseq_DmelProteins_adjusted.gff3 2> Dsim_w501_BRAKER2_RNAseq_DmelProteins_adjusted_sorting.stderr > Dsim_w501_BRAKER2_RNAseq_DmelProteins_adjusted_sorted.gff3`

The command used to generate the `Dsim_w501_BRAKER2_RNAseq_DmelProteins.gff3` input file was:

`awk 'BEGIN{FS="\t";OFS="\t";}$2=="AUGUSTUS"{if ($3 == "transcript") {split($9, genetxarr, "."); $9="gene_id \""genetxarr[1]"\"; transcript_id \""$9"\";";}; print $0;}' braker/Dsim_w501_BRAKER2/augustus.hints.gtf | ~/bin/augustus-3.3/augustus/scripts/gtf2gff.pl --printExon --gff3 --out=Dsim_w501_BRAKER2_RNAseq_DmelProteins.gff3 2> Dsim_w501_BRAKER2_RNAseq_DmelProteins_gtf2gff.stderr`

### `findSingleCopyOrthologIntrons.pl`

This Perl script uses the information provided in FASTA headers by `extractIntronsFromGFF3.pl` and a map of single-copy orthologs in order to deduce homologous introns in single-copy orthologs. Essentially, it compares the flanking exon lengths (and optionally the length of the intron itself) for a focal intron in the first species in the list, and compares these lengths pairwise to all introns in the ortholog in each remaining species, determining the index of the intron with the minimum relative difference for each length. If the indices match between the compared lengths (i.e. if the same intron has the best matches for both left and right flanking exon lengths), and both relative differences are under a threshold value (specified by `-t` as a decimal fraction), the intron is retained and output into a map TSV of a similar format to the input single-copy ortholog map.

`-f` specifies a "file of file names" (FOFN), which is just a newline-separated list of intron FASTA files, each being produced by `extractIntronsFromGFF3.pl`.

`-m` specifies the path to a single-copy ortholog map, which is a tab-separated value file consisting of the first column being an orthogroup ID, and any further columns being species-specific transcript IDs belonging to the single-copy orthogroup. The first line of this file is a header line with OG for the first column header, and species IDs for the remaining column headers. Only the content of the header before the first underscore is used to identify the intron in the input FASTA, so `extractIntronsFromGFF3.pl` should be run with the prefix (`-p`) option to match this species ID. It is perfectly acceptable to specify more than just the species ID in this prefix, `findSingleCopyOrthologIntrons.pl` should simply ignore anything after the first underscore that isn't the transcript ID.

`-t` specifies the maximum relative difference (decimal fraction, default: 0.10) allowed for flanking exon lengths to be considered a passing match.

`-i` activates the intron length filter (default is to only use flanking exon lengths)

Example Usage:

```
#Extract all the introns using only CDS features as boundaries:
extractIntronsFromGFF3.pl -i ../../refs/Dmel_ISO1.fasta -g ../../annotations/Dmel_ISO1.gff3 -p Dmel_ISO1 -b CDS | fasta_formatter -o Dmel_ISO1_all_CDSbased_introns_unwrapped.fasta
extractIntronsFromGFF3.pl -i ../../refs/Dsan_STOCAGO1482.fasta -g ../../annotations/Dsan_STOCAGO1482.gff3 -p Dsan_STOCAGO1482 -b CDS | fasta_formatter -o Dsan_STOCAGO1482_all_CDSbased_introns_unwrapped.fasta
extractIntronsFromGFF3.pl -i ../../refs/Dtei_GT53w.fasta -g ../../annotations/Dtei_GT53w.gff3 -p Dtei_GT53w -b CDS | fasta_formatter -o Dtei_GT53w_all_CDSbased_introns_unwrapped.fasta
extractIntronsFromGFF3.pl -i ../../refs/Dyak_NY73PB.fasta -g ../../annotations/Dyak_NY73PB.gff3 -p Dyak_NY73PB -b CDS | fasta_formatter -o Dyak_NY73PB_all_CDSbased_introns_unwrapped.fasta

#Prepare the FOFN of intron FASTAs:
ls *_all_CDSbased_introns_unwrapped.fasta > MSTY_CDSbased.fofn

#Identify homologous introns:
findSingleCopyOrthologIntrons.pl -f MSTY_CDSbased.fofn -m ../SCO_map.tsv -t 0.20 > MSTY_CDSbased_intron_group_map.tsv

#Partition the intron sequences into one FASTA per orthologous intron:
parseFASTARecords.pl -f MSTY_CDSbased.fofn -m MSTY_CDSbased_intron_group_map.tsv
```

### `parseFASTARecords.pl`

This Perl script partitions sequences from a set of individual-specific FASTAs into group-specific FASTAs, where the groups are groups of features, like orthologous genes, orthologous introns, etc. Only single-copy groups are currently allowed.

`-f` specifies a FOFN of FASTAs (see `findSingleCopyOrthologIntrons.pl` above), where the sequence headers in the FASTAs consist of a prefix identifying the species, then an underscore, then anything (e.g. an individual ID), then an underscore followed by the transcript ID (which should not contain underscores).

`-m` specifies a tab-separated value file mapping members of the single-copy group to a group ID (in the first column). The first line of this map file is a header specifying column names. The first column should be OG (or any other alphanumeric string, it is ignored), and the remaining column headers should be species IDs, or at least contain the species ID as a prefix before the first underscore.

The output should be a slew of FASTA files, one for each group ID in the map as long as at least one sequence matching this group was found in at least one of the input FASTAs.

For example usage, see above.

### `ScaffoldGFF3.pl`

Usage:

`ScaffoldGFF3.pl [-b <Output Broken GFF3>] -a <AGP> -g <GFF3> [-i]`

`-a` specifies an AGP file indicating the relationships between the scaffolds in the GFF3 and the chromosomes they compose.

`-g` specifies the GFF3 annotation in whichever space is input (the more-scaffolded space if `-i` is used, less-scaffolded space otherwise)

`-b` specifies the filename for the output GFF3 of features unable to be converted between coordinate spaces

`-i` is a flag meant for inverting the direction of the coordinate space transformation specified by the AGP, so if the AGP goes from scaffolds to chromosomes, using `-i` would take a chromosomal GFF3 as input, and produce a scaffold-space GFF3.

The coordinate-transformed GFF3 is output to STDOUT.

The `-d` flag can be used multiple times, and triggers output of debugging information onto STDERR.

Example Usage:

`ScaffoldGFF3.pl -dd -b Dsan_STOCAGO1482_BRAKER2_RNAseq_DmelProteins_toArms_broken.gff3 -a Dsan_STOCAGO1482_chromosome_arms_nocontam_mtDNA.agp -g Dsan_STOCAGO1482_BRAKER2_RNAseq_DmelProteins_contigs.gff3 > Dsan_STOCAGO1482_BRAKER2_RNAseq_DmelProteins_arms.gff3`

## Genome-wide statistics programs/scripts:

**Note: All C++ programs will safely compile as long as your compiler supports C++11.  I usually compile with `g++ -O3 -g -Wall --std=c++11 -o [program prefix] [program prefix].cpp`.**

### `calculatePolymorphism.cpp`

This program calculates pi given a list of FASTA filenames as positional arguments. The output columns are:

1. Scaffold ID
2. Position
3. Pi within the samples
4. Whether this site should be omitted (or, with the `-u`, is the fraction of usable sequences at that position)

Pi is calculated as in calculateDxy.

These are no longer exclusively for biallelic sites, but do reduce to the biallelic site formulae when i and j are constrained to 1 and 2.  Also, I haven't actually written out the proof yet, but I'm pretty sure that averaging the per-base values across a window is exactly equivalent to the Dxy and Pi you would calculate using Nei's formulae.

`n` in these formulae is the haploid sample size, and we assume that diploid sequences are provided (heterozygous sites identified by IUPAC degenerate bases K, M, R, S, W, and Y. Thus, without the `-i` flag, `n` is 2 times the number of input sequences.

The `-i` flag indicates that only 1 allele per site per input FASTA should be counted in the allele frequencies. The allele is chosen at random for each input sequence at each site (the `-r` argument sets the PRNG seed), and the `n` in the above formulae is equal to the number of input sequences. This is appropriate for calculations based on inbred or partially inbred individuals.

The `-s` flag causes the program to output a 1 if the site is segregating (i.e. pi > 0.0), and a 0 if not. This is the simple logical extension of `listPolyDivSites -p` to multiple samples, and can be used to evaluate Watterson's estimate of theta (you can calculate the number of segregating sites, S, sometimes denoted k, quickly from the output).

Note that you have to omit the first line (a header) in order to pass the output to `nonOverlappingWindows`, e.g. using `tail -n+2` or awk, as follows:

Further note: Versions 2.2 and up of `nonOverlappingWindows` will automatically skip the above-mentioned header line, and have the `-s` argument for selecting which column has the statistic of interest.

`calculatePolymorphism [FASTA 1] [FASTA 2] [FASTA 3] [...] | awk 'NR>1' | nonOverlappingWindows -n -w [window size in bp] -o [output TSV filename]`

or with newer versions of `nonOverlappingWindows`:

`calculatePolymorphism [FASTA 1] [FASTA 2] [FASTA 3] [...] | nonOverlappingWindows -n -w [window size in bp] -o [output TSV filename]`

In the degenerate case of inputting a single FASTA, this program behaves like `listPolyDivSites -p -n`, outputting 1 for heterozygous sites, 0 for all others.

### `decompressStats.pl`

This script adds in missing records into a subsetted statistics TSV so that the output can be used by `nonOverlappingWindows`. Any sites added in have the Omit column set to 1 so that they are not counted in the average for a window.

Usage:

`decompressStats.pl [-i stats TSV] [-f genome .fai] [-s statistic column]`

`-i` specifies the input per-site statistics TSV, with columns being Scaffold ID, Position, Statistic, and Omit. If omitted, the default is STDIN.

`-f` specifies the path to the FASTA index (.fai generated by `samtools faidx`) for the original genome.

`-s` specifies the column to use as a statistic (default is 3, but can otherwise be any integer greater than 4).

Example Usage:

`subsetVCFstats.pl -i Dyak_genomewide_Dxy.tsv -b Dyak_4fold_sites.bed | decompressStats.pl -s 3 -f Dyak.fai > Dyak_4fold_Dxy.tsv`

### `nonOverlappingWindows.cpp`

**Version change:** As of version 1.2, nonOverlappingWindows automatically skips the first line if it is a header line (i.e. does not contain any numbers), and has an option to indicate which column of the file to use as a statistic column.

Usage:

`nonOverlappingWindows [-n] [-w window size in bp] [-i input TSV] [-o output TSV] [-s stat column]`

`-n` indicates whether or not the input TSV has a fourth column indicating whether or not to omit the current site (1=omit, 0=keep).

`-w` specifies the window size to use, in base pairs.

`-i` specifies the path to the input TSV, which consists of columns: Scaffold ID, Position, Statistic, and optionally an Omit column. This argument is optional, the default is STDIN.

`-o` specifies the output TSV, which consists of columns: Scaffold ID, Position, Averaged Statistic, and an optional column indicating the fraction of used sites (i.e. 1 - fraction of omitted sites).

`-u` triggers the output of the fourth column, i.e. fraction of usable sites

`-s` specifies the column to use for averaging the statistic. The default is column 3 (i.e. `3`), but can be any integer above 2 as long as it specifies a valid column.

Example Usage:

`listPolyDivSites -n -p Dyak_NY73PB_v2_w60.fasta CY01A.fasta | nonOverlappingWindows -n -w 100000 -o CY01A_poly_w100kb.tsv`

This program takes a TSV with 3 or 4 columns, and averages the values in the 3rd column over windows of specified length.  Using the `-n` flag leads to omission of positions within the window that have a 1 in the 4th column of the input.  These sites are omitted from both the numerator and denominator of the average, hence we don't bias the estimate by interpreting masked bases as anything other than missing data.

It was originally written to calculate windowed depth using the output of `samtools depth -aa`, but the input format is general enough that most if not all of my stats tools use it.

### `subsetVCFstats.pl`

Usage:

`subsetVCFstats.pl -i [input TSV] -b [BED]`

`-i` specifies the path to the input TSV file (where the first two columns are Scaffold ID, and Position). Lines beginning with `#` will be skipped. If omitted, this argument defaults to STDIN.

`-b` specifies the BED file containing intervals to keep/subset from the input TSV.

Example Usage:

`subsetVCFstats.pl -i Dyak_Dsan_Dxy.tsv -b Dyak_4fold_sites_genomic.bed > Dyak_Dsan_Dxy_4fold_sites.tsv`

This script takes an input TSV file (where the first two columns are scaffold ID and position), and a BED file, and subsets out lines of the TSV that correspond to intervals in the BED file.  It is generally useful for subsetting, whether subsetting lines from a VCF, or Dxy or pi values from the output of `calculateDxy`, or per-base depths from the output of `samtools depth`, etc.

## MSA-related scripts:

These scripts are related to pre-processing (and post-processing) multiple sequence alignments of coding sequences. They have been used as part of a pipeline to generate MSAs of about 9,400 single-copy orthologs across Dmel, Dsan, Dtei, and Dyak, and to integrate population resequencing data into these alignments.

### `addAlignedGaps.pl`

This Perl script takes in an aligned set of CDSes from reference genomes (aligned by, say, PRANK, and provided in aligned multi-FASTA format via the `-a` argument), and an unaligned multi-FASTA of population sample CDSes (provided with the `-i` argument, or on STDIN), matches up the species prefix in their FASTA headers (e.g. `>Dsan_` in the ref alignment would match to `>Dsan_` in the population samples), and inserts gaps into the population sample sequences to align them. We perform no optimization of the MSA, simply assume that the alignment of the references is correct, and add the population samples to the output alignment (file specified by flag `-o`, or STDOUT by default). This is FAR FAR faster than inserting all the unaligned reference and population sample sequences into an MSA program, without losing much information, especially for relatively low distance alignments.

Example Usage:

`addAlignedGaps.pl -a trimmed_aligned_refs/OG0001492_prank.best.fas -i trimmed_unaligned_samples/OG0001492_trimmed.fasta -o trimmed_aligned_all/OG0001492_untrimmed.fasta`

Note: I usually make sure the input FASTAs are unwrapped, but the code should work for arbitrary line-wrapping.  I just haven't tested it thoroughly.

### `fakeHaplotype.pl`

This Perl script takes in a diploid pseudoreference FASTA (e.g. produced by the [PseudoreferencePipeline](https://github.com/YourePrettyGood/PseudoreferencePipeline/)) and randomly selects alleles at heterozygous sites to produce a haploid version of the input. The pseudo-random number generator (PRNG) seed can be set with the `-s` option (default: 42, because starting with the answer is better), and the `-a` flag allows you to choose the opposite allele, so you can run this script twice to produce the complementary two possible haploids from the input diploid. If you are splitting into both haplotypes, it is recommended to use the `-b` flag so that you can distinguish the fake haplotypes were you to concatenate the two outputs. The `-b` flag appends a `_0` or `_1` to the header of each sequence, `_0` if the `-a` flag is absent, `_1` if the `-a` flag is present. So this allows for safe splitting of haplotypes for non-haplotype-aware statistics from an outbred individual (or a synthetic diploid), and random selection of a haplotype from inbred individuals (make sure to only use one of the haplotypes, since the individual's ploidy is effectively less than 2).

Example Usage:

`cat <(fakeHaplotype.pl -i raw_CDSes_samples/OG0001492_unwrapped.fasta -b) <(fakeHaplotype.pl -i raw_CDSes_samples/OG0001492_unwrapped.fasta -a -b) | trimCDSes.awk > trimmed_unaligned_samples/OG0001492_trimmed.fasta`

This script is agnostic of line-wrapping, and will output in the same wrapping as the input.

### `skipRogers.awk`

This awk script is very very custom, and all it does is omit one of the two fake haplotypes from the Rogers et al. (2014) MBE Dyak data from the input alignment, and outputs all other sequences.

### `trimAlignment.awk`

This awk script takes in an aligned FASTA, and truncates all sequences down to the same length as the shortest sequence in the set. For the output of a typical MSA program (e.g. Clustal Omega, MUSCLE, MAFFT, PRANK), this should amount to not doing anything to the input, but this script is useful when you use an alignment of reference sequences as a guide for inserting population samples into the alignment, and the population samples happen to have differing lengths. In particular, I wrote it for use with `addAlignedGaps.pl`, where a masked pseudoreference trimmed by `trimCDSes.awk` may not truncate the terminal stop codon (because `NNN` does not match `/T(AA|AG|GA)$/`), so some of the inserted pseudoreferences have 3 extra bases. By using this script, we rectify this problem in a very very simplistic manner.

Example Usage:

`trimAlignment.awk trimmed_aligned_all/OG0001492_untrimmed.fasta > trimmed_aligned_all/OG0001492_trimmed.fasta`

### `trimCDSes.awk`

This awk script takes in an unaligned set of CDSes (for example, extracted from reference assemblies, and organized into a single FASTA file per single-copy orthogroup, as identified by OrthoFinder analysis), trims them down to length that is an integral multiple of 3, and then removes a single terminal (i.e. 3'-most) stop codon, if detected, from each sequence. At the moment, only the stop codons in the NCBI transl_table=1 genetic code are detected.

Example Usage:

`trimCDSes.awk raw_CDSes_refs/OG0001492_unwrapped.fasta > trimmed_unaligned_refs/OG0001492_trimmed.fasta`

The script is extremely short and quick, enabling extreme parallelization with GNU Parallel.

### `trimIntrons.awk`

This awk script is incredibly simple, and just trims 2 bp off of each end of an intron, to remove the constrained splice junctions.

Example Usage:

`trimIntrons.awk raw_introns_refs/OG0001492_intron1_unwrapped.fasta > trimmed_introns_refs/OG0001492_intron1_trimmed.fasta`
