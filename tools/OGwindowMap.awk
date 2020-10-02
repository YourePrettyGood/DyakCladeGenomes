#!/bin/awk -f
#Given "species" prefix, a BED of windows, and orthogroup_locations,
# this script generates a mapping from orthogroups to windows that
# contain the orthogroup.
BEGIN{
   FS="\t";
   OFS=FS;
   minoverlap=0.5;
   if (length(species) == 0) {
      print "Missing species variable, please set it." > "/dev/stderr";
      exit 1;
   };
}
#First file is a BED of windows:
FNR==NR{
   scafs[FNR]=$1;
   starts[FNR]=$2;
   ends[FNR]=$3;
   if (NF > 3) {
      names[FNR]=$4;
   };
}
#Second file is orthogroup_locations:
FNR<NR{
   #Only process orthogroup locations for the indicated species:
   split($6, proteome, "_");
   spp=proteome[1];
   if (spp == species) {
      #Iterate through all the windows for each orthogroup:
      for (i in scafs) {
         #Match on scaffold first:
         if (scafs[i] == $1) {
            #Then match on interval overlap:
            if (starts[i] < $3 && ends[i] > $2) {
               #Compute the length of the overlap:
               innerleft=starts[i] < $2 ? $2 : starts[i];
               innerright=ends[i] > $3 ? $3 : ends[i];
               overlap=innerright - innerleft + 1;
               #Compare overlap length to the gene length:
               #If the gene overlaps the window sufficiently, include it
               # in the window
               genelen=$3 - $2 + 1;
               if (overlap >= minoverlap * genelen) {
                  #Make a comma-delimited list of window IDs:
                  if (length(names[i]) > 0) {
                     id=names[i];
                  } else {
                     id=i;
                  };
                  if ($5 in windowmap) {
                     windowmap[$5]=windowmap[$5]","id;
                  } else {
                     windowmap[$5]=id;
                  };
               };
            };
         };
      };
   };
}
END{
   #Iterate through orthogroups and output window lists:
   for (i in windowmap) {
      print i, windowmap[i];
   };
}
