#!/bin/awk -f
#Given a PRNG seed and a window-orthogroup map TCSV, bootstrap each
# window and output an orthogroup-window TCSV.
#Generate the input window-orthogroup map TCSV with filterWindowOGs.awk
#Output can be fed into OGwindowAverage.awk as first file.
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(seed) == 0) {
      print "Missing PRNG seed, please set it." > "/dev/stderr";
      exit 1;
   } else {
      print "Using PRNG seed "seed > "/dev/stderr";
      srand(seed);
   };
}
#Input file is a TCSV, first column is the window ID, second column
# is a comma-delimited list of orthogroups within the window:
#Randomly draw n integers from 1 to n independently to compose
# the bootstrapped window
{
   n=split($2, ogs, ",");
   for (i=1; i<=n; i++) {
      #This is the actual bootstrap:
      ogid=int(rand()*n)+1;
      if ($1 in ogmap) {
         ogmap[$1]=ogmap[$1]","ogs[ogid];
      } else {
         ogmap[$1]=ogs[ogid];
      };
   };
}
END{
   #Convert ogmap to windowmap:
   for (i in ogmap) {
      split(ogmap[i], ogs, ",");
      for (j in ogs) {
         if (ogs[j] in windowmap) {
            windowmap[ogs[j]]=windowmap[ogs[j]]","i;
         } else {
            windowmap[ogs[j]]=i;
         };
      };
   };
   n=asorti(windowmap, ogs, "@ind_str_asc");
   #Output the orthogroups in sorted order:
   for (i=1; i<=n; i++) {
      print ogs[i], windowmap[ogs[i]];
   };
}
