#!/bin/awk -f
#Calculate assembly kmer completeness using the Merqury approach
#Arguments:
# asm: Pass-through value to identify which assembly this completeness
#      corresponds to
#The first file passed in is the output of:
# meryl statistics ${asm}.solid.meryl
# which contains the count of reliable kmers found in the assembly
#This file is computed via:
# meryl intersect output ${asm}.solid.meryl ${asm}.meryl ${reads}.gt${threshold}.meryl
#The second file passed in is the output of:
# meryl statistics ${reads}.gt${threshold}.meryl
# which contains the total count of reliable kmers
BEGIN{
   OFS="\t";
   if (length(asm) == 0) {
      print "Missing asm value, please set it." > "/dev/stderr";
      badinputs++;
   };
   if (badinputs > 0) {
      exit badinputs+1;
   };
   filenum=0;
}
FNR==1{
   filenum++;
}
filenum==1&&FNR==3{
   found=$2;
}
filenum==2&&FNR==3{
   total=$2;
}
END{
   if (length(found) == 0 || length(total) == 0) {
      print "Missing or invalid input files" > "/dev/stderr";
      exit 3;
   };
#Following Rhie et al. 2020 Genome Biology, the kmer completeness is simply
# the proportion of "reliable" kmers found in the assembly.
#"Reliable" kmers are defined as putative non-error kmers, where the
# threshold for "non-error" is determined by looking for the first
# positive slope of the kmer multiplicity histogram/density curve.
   completeness=100*found/total;
   print asm, found, total, completeness;
}
