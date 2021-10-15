#!/bin/awk -f
#This script simply extracts the OGID from the OGID tag in the fourth
# column of a BED-like file generated with SCOlocationBED.awk.
#The output is a list of OGIDs that can be used to establish a positive
# filter for all the popgen analysis.
BEGIN{
   FS="\t";
   OFS=FS;
}
{
   split($4, tags, ";");
   for (i in tags) {
      split(tags[i], elems, "=");
      if (elems[1] == "OGID") {
         print elems[2];
      };
   };
}
