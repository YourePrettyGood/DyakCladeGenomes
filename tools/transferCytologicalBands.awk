#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
}
#First file is the BED6 of SCOs for the reference with the bands tag:
FNR==NR{
   split($4, tags, ";");
   for (i in tags) {
      split(tags[i], elems, "=");
      if (elems[1] == "OGID") {
         ogid=elems[2];
      } else if (elems[1] == "bands") {
         bands=elems[2];
      };
   };
   bandlist[ogid]=bands;
}
#Second file is the BED6 of SCOs for the query assembly:
FNR<NR{
   split($4, tags, ";");
   for (i in tags) {
      split(tags[i], elems, "=");
      if (elems[1] == "OGID") {
         ogid=elems[2];
      };
   };
   print $1, $2, $3, $4";bands="bandlist[ogid], $5, $6;
}
