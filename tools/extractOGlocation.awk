#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
}
#First file is the SCO map
FNR==NR&&NR==1{
#Determine which column corresponds to the reference ID passed in:
   for (i=1;i<=NF;i++) {
      if ($i == refid) {
         refcol=i;
      }
   }
}
FNR==NR&&NR>1{
#Read in the orthogroup mappings for the specified reference:
   SCOmap[$refcol]=$1;
}
#Second file is the annotation for the specified reference
FNR<NR&&$3=="mRNA"{
#Output coordinates of the transcript if it's found in the SCO map:
   split($9, tags, ";");
   for (tag in tags) {
      split(tags[tag], elems, "=");
      if (elems[1] == "ID") {
         if (elems[2] in SCOmap) {
            print $1, $4, $5, $7, SCOmap[elems[2]], refid, elems[2];
         }
      }
   }
}
