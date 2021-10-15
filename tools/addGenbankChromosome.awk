#!/bin/awk -f
#This script adds the Genbank chromosome specifier to the
# scaffolds of the input FASTA, excluding those with IDs
# matching a provided pattern.
#If an mtDNA pattern is provided, the appropriate specifiers
# for the mtDNA contig are added, and the mtDNA name is
# simplified to mtDNA to avoid name length errors
BEGIN{
   #Blank pattern or no pattern means to augment all headers
}
#Add specifiers to the header lines:
/^>/{
   shortid=substr($1, 2);
   if (length(exclude) > 0 && $0 ~ exclude) {
      print $0" [gcode=1]";
   } else if (length(mtDNA) > 0 && $0 ~ mtDNA) {
      print ">mtDNA [chromosome=MT] [location=mitochondrion] [topology=circular] [completeness=complete] [mgcode=5]"
   } else {
      print $0" [chromosome="shortid"] [gcode=1]";
   };
}
#Feed through sequence lines:
!/^>/{
   print;
}
