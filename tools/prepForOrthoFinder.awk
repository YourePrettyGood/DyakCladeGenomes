#!/bin/awk -f
#Rename proteins to work as input for OrthoFinder, and make a map
# to convert IDs back.
BEGIN{
   i=0;
   mapfn=species".map";
   speciespattern=">"species"_";
}
#Just feed through any sequence lines:
!/^>/{
   print;
}
/^>/{
#Rename protein to [species ID]_#
   print ">"species"_"i;
   txid=$0;
#Trim the species prefix from the map:
   sub(speciespattern, "", txid);
#Remove the translation frame appended by transeq:
   sub(/_1$/, "", txid);
#Output to the map file:
   print txid"\t"species"_"i > mapfn;
   i+=1;
}
