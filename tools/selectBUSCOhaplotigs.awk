#!/bin/awk -f
#This script takes in the STDOUT output of scoreBUSCOhaplotigs.awk
# and applies a set of rules to determine which contigs/scaffolds
# to purge/select.
#The rules are:
#1) Column 4 (DupRateA) OR Column 8 (DupRateB) must be 1.0
#   i.e. all BUSCO hits on scaffold A OR scaffold B must be of
#   type "Duplicated" (the OR is not exclusive)
#2) Column 10 (SharedDupRateA) OR Column 11 (SharedDupRateB) must
#   be 1.0
#   i.e. all BUSCO hits on scaffold A OR scaffold B must be shared
#   with the other scaffold in the pair
#3a)Given 1 and 2 pass, if Col 4 == Col 8 AND Col 10 == Col 11,
#   pick the shorter scaffold (by length) for purging
#   (in a tie, pick Scaffold A)
#3b)Given 1 and 2 pass, pick the scaffold with DupRate == 1.0 for
#   purging
#Optionally, the corresponding line of the input can be appended
# to the purged scaffold ID for debugging when the debug option
# is passed
#
#Note about output:
#The output may contain redundant scaffold IDs if the relationships
# are more complicated than simply primary-alternate. In testing,
# I saw a pair of alternate haplotigs of different lengths that had
# exactly the same number of BUSCO hits, and were both alternate to
# a major scaffold, thus they both need to be purged. The shorter
# alternate was listed twice. Thus, check the output with debug,
# but you may be able to just take the non-debug output and pipe
# through sort | uniq.
BEGIN{
   #Set input and output delimiters to tabs:
   FS="\t";
   OFS=FS;
}
#Input file columns:
#1) Scaffold A ID
#2) Scaffold A Duplicated BUSCO count
#3) Scaffold A Total BUSCO count
#4) Scaffold A Duplicated BUSCO rate ($2/$3)
#5) Scaffold B ID
#6) Scaffold B Duplicated BUSCO count
#7) Scaffold B Total BUSCO count
#8) Scaffold B Duplicated BUSCO rate ($6/$7)
#9) Shared BUSCO count (necessarily only Duplicated BUSCOs)
#10)Shared BUSCO rate A ($9/$3)
#11)Shared BUSCO rate B ($9/$7)
#12)Scaffold A length
#13)Scaffold B length
#Skip the header line and check for Rule 1:
NR>1&&$4==1.0||$8==1.0{
   #Apply Rule 2:
   if ($10 == 1.0 || $11 == 1.0) {
      #Rule 3a:
      if ($4 == $8 && $10 == $11) {
         if ($12 > $13) { #Using > here ensures outputting Scaffold A in a tie
            if (length(debug) > 0) {
               print $5, $0;
            } else {
               print $5;
            };
         } else {
            if (length(debug) > 0) {
               print $1, $0;
            } else {
               print $1;
            };
         };
      #Rule 3b:
      } else {
         if ($4 == 1.0) { #We only need to check one because Rule 1 is OR
            if (length(debug) > 0) {
               print $1, $0;
            } else {
               print $1;
            };
         } else {
            if (length(debug) > 0) {
               print $5, $0;
            } else {
               print $5;
            };
         };
      };
   };
}
