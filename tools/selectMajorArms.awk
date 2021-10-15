#!/bin/awk -f
#A simple script to only select the major chromosome arms from a BED file
BEGIN{
   FS="\t";
   OFS=FS;
}
{
#   if ($1 ~ /^(Scf_)?[X234][LR]?/) {
   if ($1 ~ /^(Scf_)?[X23][LR]?/) {
      print;
   };
}
