#!/bin/awk -f
#chromlist should be specified as a comma-separated list of chromosomes
# in the order they should appear
#Parse the chromlist and put it into a hash:
BEGIN{
   FS="\t";
   OFS=FS;
   split(chromlist, chroms, ",");
   for (i in chroms) {
      chromhash[chroms[i]]=1;
   }
}
#Extract the length of the chromosomes in chromlist from the input FAI:
$1 in chromhash{
   chromlen[$1]=$2;
}
#Output the plotting order in the order of the input list:
END{
   for (chrom in chroms) {
      if (chromlen[chroms[chrom]]>0) {
#We force everything in the same orientation here
         print chroms[chrom], chromlen[chroms[chrom]], "+";
      }
   }
}
