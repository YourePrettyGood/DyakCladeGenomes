#!/bin/awk -f
#Note: an "id" variable will append an extra column to the output.
#Given "species" prefix, statistic column name, weight column name,
# an orthogroup-window map TCSV file, and filtered Polymorphorama
# output file with header line, this script averages the statistic
# across orthogroups within each window, performing a weighted
# average if the weightcol variable matches a column in the header.
#Use OGwindowMap.awk to generate the orthogroup-window map TCSV.
#The filtered Polymorphorama output file's header must include
# two particular column names: OG and Species
#These column names are case-sensitive, and represent the orthogroup
# and species for the statistics on that line. OG is NOT outgroup.
BEGIN{
   FS="\t";
   OFS=FS;
   #weight of 0 means take a naive average:
   weight=0;
   if (length(statcol) == 0) {
      print "Missing statistic column name, please set it." > "/dev/stderr";
      exit 1;
   }
   if (length(species) == 0) {
      print "Missing species variable, please set it." > "/dev/stderr";
      exit 2;
   };
   if (length(id) == 0) {
      print "id column of output will be blank." > "/dev/stderr";
   } else {
      print "id column will show "id > "/dev/stderr";
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
#The weightcol variable may match a numerical column's name, or be omitted
FNR<NR&&FNR==1{
   #Find column indexes by matching column names in header line:
   for (i=1; i<=NF; i++) {
      if ($i == statcol) {
         stat=i;
      };
      if ($i == weightcol) {
         weight=i;
      };
      if ($i == "OG") {
         og=i;
      };
      if ($i == "Species") {
         spp=i;
      };
   };
   if (length(debug) > 0) {
      print "OG="og",Species="spp","statcol"="stat","weightcol"="weight > "/dev/stderr";
   };
}
FNR<NR&&FNR>1&&$spp==species{
   split(windowmap[$og], windowids, ",");
   for (i in windowids) {
      #Omit any orthogroups with NA statistic values:
      if ($stat != "NA") {
         #Keep track of the windows (we'll sort them at the end):
         windows[windowids[i]]++;
         #Default is to use weight of 1 for naive average:
         w=1;
         #If the weight column was set and found, use weighted average:
         if (weight > 0) {
            if (length(debug) > 0) {
               print $og, windowids[i], weight"="$weight > "/dev/stderr";
            };
            w=$weight;
         };
         #Weighted average is done via sum(x_i*w_i)/sum(w_i)
         numerator[windowids[i]]+=$stat * w;
         denominator[windowids[i]]+=w;
         if (length(debug) > 0) {
            print $og, windowids[i], $stat, $stat*w, w, id > "/dev/stderr";
         };
      };
   };
}
END{
   n=asorti(windows, sortedwindows, "@ind_num_asc");
   if (length(noheader)==0) {
      #Output a header line:
      print "WindowID", "Numerator", "Denominator", "Average"statcol, "NumOGs", "ID";
   }
   #Output the window averages in sorted order:
   for (i=1; i<=n; i++) {
      if (sortedwindows[i] in numerator && sortedwindows[i] in denominator) {
         if (denominator[sortedwindows[i]] == 0) {
            print sortedwindows[i], numerator[sortedwindows[i]], denominator[sortedwindows[i]], "NA", windows[sortedwindows[i]], id;
         } else {
            print sortedwindows[i], numerator[sortedwindows[i]], denominator[sortedwindows[i]], numerator[sortedwindows[i]]/denominator[sortedwindows[i]], windows[sortedwindows[i]], id;
         };
      } else {
         print sortedwindows[i], "NA", "NA", "NA", "NA", id;
      };
   };
}
