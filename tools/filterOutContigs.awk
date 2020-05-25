#!/bin/awk -f
#For the time being, we make the assumption that the contig IDs in the
# first file match the first substring after the > and before a space
# in the FASTA
#The script's default function is to filter contigs from a FASTA, but
# given -v "inputtype=GFF" or -v "inputtype=GTF" (case-insensitive)
BEGIN{
   keep=1;
#Default input is FASTA:
   if (length(inputtype) == 0) {
      inputtype="FASTA";
   }
#Make sure to parse and output tab-separated if GTF or GFF is input:
   if (tolower(inputtype) ~ /g[ft]f3?/) {
      FS="\t";
      OFS=FS;
   }
}
#First file is the list of contigs to filter out:
FNR==NR{
   filterset[$1]=1;
}
#Second file is the FASTA of contigs (or a GFF or GTF):
FNR<NR{
   if (tolower(inputtype) ~ /fa(sta)?/) {
#Condition on header line not matching filterset:
      if ($0 ~ /^>/) {
         keep=1;
         contigname=substr($1, 2);
         if (length(filterset[contigname]) > 0) {
            keep=0;
         }
      }
#General case for GFF or GTF:
   } else if (tolower(inputtype) ~ /g[ft]f3?/) {
      keep=1;
#Special case of GFF3, we need to eliminate ##sequence-region headers
# for filtered contigs:
      if (tolower(inputtype) ~ /gff3?/) {
         if ($0 ~ /^##sequence-region/) {
            split($0, headerparts, " ");
            if (length(filterset[headerparts[2]]) > 0) {
               keep=0;
            }
         }
      }
#Otherwise, filter non-header records:
      if (!($0 ~ /^#/)) {
         if (length(filterset[$1]) > 0) {
            keep=0;
         }
      }
   } else {
      print "Unknown inputtype "inputtype", exiting" > "/dev/stderr";
      exit 1;
   }
#Print if not filtered out:
   if (keep) {
      print;
   }
}
