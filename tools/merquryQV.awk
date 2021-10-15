#!/bin/awk -f
#Calculate assembly QV using the Merqury approach
#Arguments:
# k:   kmer length, which is a necessary component of the error calculation
# asm: Pass-through value to identify which assembly this QV corresponds to
#The first file passed in is the output of:
# meryl statistics ${asm}.0.meryl
# which contains the count of kmers found in the assembly but not the reads
#The second file passed in is the output of:
# meryl statistics ${asm}.meryl
# which contains the total count of kmers found in the assembly
BEGIN{
   OFS="\t";
   if (length(k) == 0) {
      print "Missing k value, please set it." > "/dev/stderr";
      badinputs++;
   };
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
filenum==1&&FNR==4{
   asmonly=$2;
}
filenum==2&&FNR==4{
   total=$2;
}
END{
   if (length(asmonly) == 0 || length(total) == 0) {
      print "Missing or invalid input files" > "/dev/stderr";
      exit 4;
   };
#Following Rhie et al. 2020 Genome Biology's adaptation of the Ondov et al.
# 2019 Genome Biology containment score, the error probability is estimated
# using a binomial kmer survival model, assuming that the read population
# contains all kmers in the assembly. In that case, the probability of
# observing all k kmers that cover a given base is (p_match)^k, which is
# approximated by K_shared/K_total when averaging across bases and assuming
# kmers are unique and independent. Thus we have the relation:
# (p_match)^k = K_shared/K_total
#Rearranging, we have:
# \hat{p_match} = (K_shared/K_total)^(1/k)
#Noting that K_total = K_shared + K_asmonly and that we observe K_asmonly:
# \hat{p_match} = (1-K_asmonly/K_total)^(1/k)
#Now this is a match probability across all bases, so we can estimate
# an error probability by simply taking the additive inverse:
# \hat{p_error} = 1-(1-K_asmonly/K_total)^(1/k) 
   error=1-(1-asmonly/total)^(1/k);
#For the QV, we just PHRED-scale the error probability, accounting for
# the fact that awk's log() function computes the natural log, so we
# use the change of base formula:
   qv=-10*log(error)/log(10);
   print asm, asmonly, total, qv, error;
}
