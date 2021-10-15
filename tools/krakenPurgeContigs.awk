#!/bin/awk -f
#This is a very simple script to output the contigs/scaffolds from the
# STDOUT of a kraken2 run with --use-names, identifying any contigs
# that have consensus taxonomic classification other than Drosophila or
# Homo. The intended kraken2 database for this is composed of the
# "Standard" database (Human, RefSeq archaea, bacteria, viruses,
# UniVec), and then added the RefSeq assemblies of Drosophila melanogaster
# and D. yakuba.
BEGIN{
   FS="\t";
   OFS=FS;
}
$1=="C"&&$3!~/^(Drosophila|Homo)/{
   print $2;
}
