#!/bin/awk -f
#Assumptions:
#Input BUSCO full_table*.tsv is in sorted BUSCO id order, and
# hits to major chromosomes/arms come first before haplotig hits.
BEGIN{
   #Parse and output tab-separated:
   FS="\t";
   OFS=FS;
   #Default for major arms is X,2L,2R,3L,3R,4:
   if (length(major) == 0) {
      major="X,2L,2R,3L,3R,4";
   }
   #Split the comma-separated string of major chromosomes/arms:
   split(major, a, ",");
   #Make them keys in an array for fast existence checks:
   for (i in a) {
      arms[a[i]]=1;
   };
   #Keep track of previous BUSCO id to output on transition:
   previd="";
   #Make sure array iteration order is consistent:
   PROCINFO["sorted_in"]="@val_str_asc";
   #Other options include:
   #merge: Change the primary hit to be Complete instead of Duplicated
   #purge_nonmajor: Output secondary hits for BUSCOs where the primary hit
   # is to a nonmajor scaffold
   #Thoughts to perhaps add later:
   #thresh: Maximum scaffold length relative to primary hit
}
#For thresh, would need this:
##First file is the FAI or .genome for scaffold lengths:
#FNR==NR{
#   scaflen[$1]=$2;
#}
##Second file is the full_table_*.tsv from BUSCO in the run_* directory
#File is the full_table_*.tsv from BUSCO in the run_* directory
#Feed through commented lines from the second file to STDERR:
#FNR<NR&&/^#/{
/^#/{
   print > "/dev/stderr";
}
#For non-commented lines, screen for Duplicated BUSCOs and deduplicate:
#FNR<NR&&!/^#/{
!/^#/{
   #If previous BUSCO id is different from current, output any
   # duplicated (or deduplicated) records:
   if ($1 != previd) {
      for (i in duprecs) {
         #If only 1 record remains, convert to Complete if requested:
         if (length(duprecs[2]) == 0 && length(merge) > 0) {
            sub("Duplicated", "Complete", duprecs[i]);
         };
         print duprecs[i] > "/dev/stderr";
      };
      delete duprecs;
   };
   #Feed through any non-Duplicated BUSCOs to STDERR:
   if ($2 != "Duplicated") {
      print > "/dev/stderr";
   } else {
      #Take care of purging nonmajor-nonmajor duplicates if requested:
      if (length(purge_nonmajor) > 0) {
         nonmajor = 1;
      } else {
         nonmajor = dup[$1] in arms;
      };
      #Keep the first hit if it's on a major chromosome/arm:
      if ($3 in arms) {
         #Store the major copy's scaffold:
         dup[$1]=$3;
         #Store the record unmodified for now:
         duprecs[length(duprecs)+1]=$0;
#         if (length(merge) > 0) {
#            #Merge by retaining this hit and changing from Duplicated to Complete:
#            $2="Complete";
#         }
#         print > "/dev/stderr";
      } else if (length(dup[$1]) > 0 && nonmajor) {
         #Make sure not to output multiple copies of the same scaffold:
         if (length(nodup[$3]) == 0) {
            #Print the scaffold to STDOUT if it's a major-nonmajor duplicate:
            print $3;
         };
         nodup[$3]=1;
      } else {
         #Store the first copy's scaffold:
         dup[$1]=$3;
#         #Feed the line through to STDERR since this is a nonmajor-nonmajor
#         # duplicate:
#         print > "/dev/stderr";
         duprecs[length(duprecs)+1]=$0;
      };
   };
   previd=$1;
}
END{
   for (i in duprecs) {
      #If only 1 record remains, convert to Complete if requested:
      if (length(duprecs[2]) == 0 && length(merge) > 0) {
         sub("Duplicated", "Complete", duprecs[i]);
      };
      print duprecs[i] > "/dev/stderr";
   };
   delete duprecs;
}
