#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(excludedcols)==0) {
      excludedcols="";
   };
   split(excludedcols, excludedcolarr, " ");
   split("", rminds);
}
NR==1{
#Only keep parts of the headers before the first underscore:
   for (i=2; i<=NF; i++) {
#Identify the columns to be skipped:
      for (j in excludedcolarr) {
         if ($i ~ excludedcolarr[j]) {
            rminds[j]=i;
         };
      };
#Trim off rest of headers after underscore:
      split($i, fieldarr, "_");
      $i=fieldarr[1];
   };
   asort(rminds, sortedrminds, "@val_num_desc");
#Remove from back to front to avoid index shifts with each removal:
   for (i in sortedrminds) {
      $0=gensub(/\s*\S+/, "", sortedrminds[i]);
   };
   print $0;
}
NR>1{
#Skip columns again from back to front:
   for (i in sortedrminds) {
      $0=gensub(/\s*\S+/, "", sortedrminds[i]);
   };
   print $0;
}
