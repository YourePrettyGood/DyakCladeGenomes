#!/bin/awk -f
#chromlist should be specified as a comma-separated list of chromosomes
# in the order they should appear
#If supplied, nonchrom indicates that only unscaffolded contigs should
# be included in the output
#Parse the chromlist and put it into a hash:
BEGIN{
   FS="\t";
   OFS=FS;
   split(chromlist, chroms, ",");
   for (i in chroms) {
      chromhash[chroms[i]]=1;
   }
   if (length(nonchrom) == 0) {
      nonchrom=0;
   } else {
      nonchrom=1;
   }
}
#First file passed in is the scaffolding config string, so parse it
FNR==NR{
   split($0, scafs, "=>");
   for (i=1; i<=length(scafs); i++) {
      split(scafs[i], ctgs, "->");
      split(ctgs[1], scaf, ":");
      scafname=scaf[1];
      ctgs[1]=scaf[2];
      ctgnum[scafname]=length(ctgs);
      for (j=1; j<=ctgnum[scafname]; j++) {
         if (substr(ctgs[j], length(ctgs[j])) == "*") {
            orientation="-";
            ctgs[j]=substr(ctgs[j], 1, length(ctgs[j])-1);
         } else {
            orientation="+";
         }
         scaffoldedctgs[ctgs[j]]=1;
         ctgorder[scafname, j]=ctgs[j];
         ctgorient[scafname, j]=orientation;
      }
   }
}
#Extract the lengths of the contigs from the input FAI (the second file):
FNR<NR{
   ctglen[$1]=$2;
   if (!($1 in scaffoldedctgs)) {
      nonchromctgs[$1]=1;
   }
}
#Output the plotting order in the order of the input list:
END{
   if (nonchrom) {
      PROCINFO["sorted_in"]="@ind_str_asc";
      for (i in nonchromctgs) {
         print i, ctglen[i], "+";
      }
   } else {
      for (i in chroms) {
         chrom=chroms[i];
         if (ctgnum[chrom]>0) {
            for (j=1; j<=ctgnum[chrom]; j++) {
               ctg=ctgorder[chrom, j];
               print ctg, ctglen[ctg], ctgorient[chrom, j];
            }
         }
      }
   }
}
