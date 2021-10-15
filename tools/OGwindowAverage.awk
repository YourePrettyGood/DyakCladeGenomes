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
#Revisions on 2020/08/30 enable passing a comma-separated list of
# statistic column names. If a list is passed, then the output header
# names change somewhat:
# Numerator becomes Sum[statname], Denominator becomes SumWeight,
# and there is an Average[statname] column for each statistic.
#These appear in order, iterating over statistics, and NumOGs and
# ID remain the last column names.
#Revision on 2021/07/12 allows for indicating the species and OG
# column names with the sppcol and ogcol variables
# (backwards compatibility via defaults of Species and OG)
BEGIN{
   FS="\t";
   OFS=FS;
   #weight of 0 means take a naive average:
   weight=0;
   if (length(statcol) == 0) {
      print "Missing statistic column name, please set it." > "/dev/stderr";
      exit 1;
   }
   n_stats=split(statcol, statcols, ",");
   if (length(species) == 0) {
      print "Missing species variable, please set it." > "/dev/stderr";
      exit 2;
   };
   if (length(id) == 0) {
      print "id column of output will be blank." > "/dev/stderr";
   } else {
      print "id column will show "id > "/dev/stderr";
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
#The weightcol variable may match a numerical column's name, or be omitted
FNR<NR&&FNR==1{
   #Find column indexes by matching column names in header line:
   for (i=1; i<=NF; i++) {
      for (j=1; j<=n_stats; j++) {
         if ($i == statcols[j]) {
            stats[j]=i;
         };
      };
      if ($i == weightcol) {
         weight=i;
      };
      if ($i == ogcol) {
         og=i;
      };
      if ($i == sppcol) {
         spp=i;
      };
   };
   if (length(debug) > 0) {
      printf "OG=%s,Species=%s", og, spp > "/dev/stderr";
      for (j=1; j<=n_stats; j++) {
         printf ",%s=%.6g", statcols[j], stats[j] > "/dev/stderr";
      };
      printf ",%s=%.6g\n", weightcol, weight > "/dev/stderr";
#      print "OG="og",Species="spp","statcol"="stat","weightcol"="weight > "/dev/stderr";
   };
}
FNR<NR&&FNR>1&&$spp==species{
   split(windowmap[$og], windowids, ",");
   for (i in windowids) {
      #Default is to use weight of 1 for naive average:
      w=1;
      #If the weight column was set and found, use weighted average:
      if (weight > 0) {
         if (length(debug) > 0) {
            print $og, windowids[i], weight"="$weight > "/dev/stderr";
         };
         w=$weight;
      };
      addwindow=0;
      for (j=1; j<=n_stats; j++) {
         #Omit any orthogroups with NA statistic values:
         if ($stats[j] != "NA") {
            #Keep track of the windows (we'll sort them at the end):
            addwindow=1;
            #Weighted average is done via sum(x_i*w_i)/sum(w_i)
            numerator[windowids[i],j]+=$stats[j] * w;
            denominator[windowids[i],j]+=w;
            if (length(debug) > 0) {
               print $og, windowids[i], $stats[j], $stats[j]*w, w, id > "/dev/stderr";
            };
         };
      };
      #Keep track of the windows (we'll sort them at the end):
      if (addwindow) {
         windows[windowids[i]]++;
      };
   };
}
END{
   n=asorti(windows, sortedwindows, "@ind_num_asc");
   if (length(noheader)==0) {
      #Output a header line:
      #Special header if only 1 stat column:
      if (n_stats == 1) {
         print "WindowID", "Numerator", "Denominator", "Average"statcols[1], "NumOGs", "ID";
      } else {
         printf "WindowID"OFS;
         for (j=1; j<=n_stats; j++) {
            printf OFS"Sum%s"OFS"SumWeight%s"OFS"Average%s", statcols[j], statcols[j], statcols[j];
         };
         printf OFS"NumOGs"OFS"ID\n";
      };
   }
   #Output the window averages in sorted order:
   for (i=1; i<=n; i++) {
      printf "%s", sortedwindows[i];
      for (j=1; j<=n_stats; j++) {
         if ((sortedwindows[i],j) in numerator && (sortedwindows[i],j) in denominator) {
            if (denominator[sortedwindows[i],j] == 0) {
               printf OFS"%.6g"OFS"%.6g"OFS"NA", numerator[sortedwindows[i],j], denominator[sortedwindows[i],j];
#               print sortedwindows[i], numerator[sortedwindows[i]], denominator[sortedwindows[i]], "NA", windows[sortedwindows[i]], id;
            } else {
               printf OFS"%.6g"OFS"%.6g"OFS"%.6g", numerator[sortedwindows[i],j], denominator[sortedwindows[i],j], numerator[sortedwindows[i],j]/denominator[sortedwindows[i],j];
#               print sortedwindows[i], numerator[sortedwindows[i]], denominator[sortedwindows[i]], numerator[sortedwindows[i]]/denominator[sortedwindows[i]], windows[sortedwindows[i]], id;
            };
         } else {
            printf OFS"NA"OFS"NA"OFS"NA";
         };
      };
      printf OFS"%u"OFS"%s\n", windows[sortedwindows[i]], id;
   };
}
