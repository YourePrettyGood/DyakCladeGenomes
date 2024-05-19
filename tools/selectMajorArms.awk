#!/bin/awk -f
#A simple script to only select the major chromosome arms from a BED file
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(inputtype) == 0) {
      inputtype="BED";
   };
   print "Input assumed to be "inputtype > "/dev/stderr";
   #Define the regex for major chromosome arms here:
   #majorarms="^(Scf_)?[X234][LR]?";
   majorarms="^(Scf_)?[X23][LR]?";
   #Set the flag for keeping records to skip by default:
   keep=0;
}
inputtype=="BED"{
   if ($1 ~ majorarms) {
      print;
   };
}
inputtype=="FASTA"&&/^>/{
   if (substr($1, 2) ~ majorarms) {
      keep=1;
      print;
   } else {
      keep=0;
   };
}
inputtype=="FASTA"&&!/^>/&&keep==1{
   print;
}
