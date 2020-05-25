#!/bin/awk -f
#Generates a gene_transcript_map for TransDecoder from a GFF3
# based on the Parent field of mRNA records
BEGIN{
#Make sure we parse and output tab-separated:
   FS="\t";
   OFS="\t";
}
NF>=9 && $3 == "mRNA"{
   split($9, tagarr, ";");
   for (tag in tagarr) {
      split(tagarr[tag], elems, "=");
      if (elems[1] == "ID") {
         tx=elems[2];
      } else if (elems[1] == "Parent") {
         gene=elems[2];
      };
   };
   print gene, tx;
}
