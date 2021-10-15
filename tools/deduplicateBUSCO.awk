#!/bin/awk -f
#This script deduplicates the BUSCO full table for re-analysis
# after purging duplicate and contaminated scaffolds. The first
# input file is a TCSV consisting of a column of scaffold names
# and a second column being a comma-separated list of BUSCOs
# found on that scaffold. The TCSV is filtered to only retain
# those scaffolds that survive the purge. The second input file
# is the original BUSCO full table which will be filtered.
BEGIN{
   #Parse and output tab-separated:
   FS="\t";
   OFS=FS;
   #Make sure array iteration order is consistent:
   PROCINFO["sorted_in"]="@ind_str_asc";
}
#First file is a TCSV with two columns:
# 1) Scaffold name
# 2) Comma-separated list of BUSCOs on the scaffold
#Any purged scaffolds will be absent from this file.
FNR==NR{
   split($2, buscos, ",");
   for (b in buscos) {
      if (length(busco_loc[buscos[b]]) == 0) {
         busco_loc[buscos[b]]=$1;
      } else {
         busco_loc[buscos[b]]=busco_loc[buscos[b]]","$1;
      };
   };
}
##Second file is the full_table_*.tsv from BUSCO in the run_* directory
#File is the full_table_*.tsv from BUSCO in the run_* directory
#Feed through commented lines from the second file:
FNR<NR&&/^#/{
#/^#/{
   print;
}
#For non-commented lines, screen for Duplicated BUSCOs and deduplicate:
FNR<NR&&!/^#/{
#!/^#/{
   #First check if Missing and feed through if so:
   if ($2 == "Missing") {
      print;
   } else if ($1 in busco_loc) { #Skip any BUSCOs not found in the filtered list
      n=split(busco_loc[$1], a, ",");
      if ($2 == "Duplicated" && n == 1) { #Adjust the category for a deduplicated BUSCO
         #print "Deduplicated "$1" on scaffold "$3 > "/dev/stderr";
         $2="Complete";
      };
      #Only print the (adjusted) line if the location is in the purged list:
      for (i in a) {
         if ($3 == a[i]) {
            print $0;
         };
      };
   } else {
      #print $1" not found in busco_loc hash" > "/dev/stderr";
   };
}
