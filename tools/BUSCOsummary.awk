#!/bin/awk -f
#This script produces a summary of the BUSCO results given the
# full_table_*.tsv file in the BUSCO run_* directory.
#The default is to produce the standard BUSCO text results string, but
# by specifying the Rout flag, the results will be output as a tall
# table in tab-separated format.
#Specifying the "species" variable will print an extra column with that value.
#This option is useful for combining results across multiple assemblies.
#Added an option 2021/01/24 to output in a more R-readable format if
# requested
#Without this option, the default is to output BUSCO-formatted strings
#Also added an option to prepend the species name in Rout mode
BEGIN{
   #Read in and output TSV:
   FS="\t";
   OFS=FS;
   #Make sure to initialize counts to 0:
   counts["S"]=0;
   counts["D"]=0;
   counts["F"]=0;
   counts["M"]=0;
}
!/^#/{
   #Partition by BUSCO type:
   if ($2 == "Duplicated") {
      #Make sure not to double-count Duplicated BUSCOs:
      if ($1 in nodup) {
      } else {
         counts["D"]+=1;
      };
      nodup[$1]=1;
   } else if ($2 == "Complete") {
      counts["S"]+=1;
   } else if ($2 == "Fragmented") {
      counts["F"]+=1;
   } else if ($2 == "Missing") {
      counts["M"]+=1;
   } else {
      print "Unknown Type", $0;
   };
}
END{
   #Complete = Single + Duplicated:
   counts["C"]=counts["S"]+counts["D"];
   #Total = Complete + Fragmented + Missing:
   counts["N"]=counts["C"]+counts["F"]+counts["M"];
   #Rout format is tall: Category, Count, Percentage
   if (length(noheader) == 0) {
      #Print a header line:
      if (length(species) > 0) {
         printf "Species%s", OFS;
      };
      print "Category", "Count", "Percentage";
      #Print out the counts and percentages (order is deterministic, lexicographical):
      PROCINFO["sorted_in"]="@ind_str_asc";
      for (c in counts) {
         if (length(species) > 0) {
            printf "%s%s", species, OFS;
         };
         print c, counts[c], counts[c]*100/counts["N"];
      };
   } else {
      #Print the counts:
      print "C:"counts["C"]"[S:"counts["S"]",D:"counts["D"]"],F:"counts["F"]",M:"counts["M"]",n="counts["N"];
      #Print the percentages:
      print "C:"counts["C"]*100/counts["N"]"%[S:"counts["S"]*100/counts["N"]"%,D:"counts["D"]*100/counts["N"]"%],F:"counts["F"]*100/counts["N"]"%,M:"counts["M"]*100/counts["N"]"%,n="counts["N"];
   };
}
