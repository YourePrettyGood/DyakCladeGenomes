#!/bin/awk -f
#Given "species" prefix, statistic column name, an orthogroup-window
# map TCSV file, and filtered Polymorphorama output file with header
# line, this script filters the orthogroup-window map and converts
# it to its inverse map. This inverse map facilitates bootstrapping.
#Use OGwindowMap.awk to generate the orthogroup-window map TCSV.
#The filtered Polymorphorama output file's header must include
# two particular column names: OG and Species
#These column names are case-sensitive, and represent the orthogroup
# and species for the statistics on that line. OG is NOT outgroup.
#Revisions on 2021/07/12 to specify the names of the orthogroup and
# species header names via the ogcol and sppcol variables
# (backwards compatibility maintained by setting defaults of Species and OG)
BEGIN{
   FS="\t";
   OFS=FS;
   if (length(statcol) == 0) {
      print "Missing statistic column name, please set it." > "/dev/stderr";
      exit 1;
   }
   if (length(species) == 0) {
      print "Missing species variable, please set it." > "/dev/stderr";
      exit 2;
   };
   if (length(sppcol) == 0) {
      sppcol="Species";
   };
   if (length(ogcol) == 0) {
      ogcol="OG";
   };
}
#First file is a TCSV, first column is the orthogroup, second column
# is a comma-delimited list of window IDs containing the orthogroup:
FNR==NR{
   windowmap[$1]=$2;
}
#Second file is a filtered version of the statistics output by
# Polymorphorama (not direct output):
#Two required columns are "OG" and "Species", indicating the orthogroup
# and ingroup species represented on the current line.
#The value of the statcol variable must match a numerical column's name
FNR<NR&&FNR==1{
   #Find column indexes by matching column names in header line:
   for (i=1; i<=NF; i++) {
      if ($i == statcol) {
         stat=i;
      };
      if ($i == ogcol) {
         og=i;
      };
      if ($i == sppcol) {
         spp=i;
      };
   };
   if (length(debug) > 0) {
      print "OG="og",Species="spp","statcol"="stat > "/dev/stderr";
   };
}
FNR<NR&&FNR>1&&$spp==species{
   split(windowmap[$og], windowids, ",");
   for (i in windowids) {
      #Omit any orthogroups with NA statistic values:
      if ($stat != "NA") {
         if (windowids[i] in ogmap) {
            ogmap[windowids[i]]=ogmap[windowids[i]]","$og;
         } else {
            ogmap[windowids[i]]=$og;
         };
         if (length(debug) > 0) {
            print $og, windowids[i] > "/dev/stderr";
         };
      };
   };
}
END{
   n=asorti(ogmap, windows, "@ind_num_asc");
   #Output the windows in sorted order:
   for (i=1; i<=n; i++) {
      print windows[i], ogmap[windows[i]];
   };
}
