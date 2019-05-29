#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
}
NR==1{
#Adjust header line so column 1 header is "OG", and store number of species
# in header:
   $1="OG";
   print $0;
   num_spp=NF;
}
NR>1{
#We use a simple trick to identify single-copy orthogroups:
# If the product of the number of genes in the orthogroup across species
# is 1, then each element of the product must have been 1, thus SCO.
   sco=1;
   for (i=2; i<=NF; i++) {
      sco=sco*split($i, txarr, ",");
   };
#We also have to check that all species are represented in the SCO,
# hence NF==num_spp:
   if (sco==1 && NF==num_spp) {
      print;
   };
}
