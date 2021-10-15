#!/bin/awk -f
#Assumptions:
#Input BUSCO full_table*.tsv is in sorted BUSCO id order
#This script scores contigs/scaffolds based on their quantity
# and fraction of Duplicated BUSCOs, and how many Duplicated
# BUSCOs they share with other contigs/scaffolds. These scores
# can be used to decide which contigs/scaffolds to filter out,
# and thus which BUSCOs are deduplicated.
#It also outputs a comma-separated list of BUSCO IDs for each
# scaffold on STDERR so that we can use this for deduplication
# alongside a list of contigs/scaffolds to retain.
BEGIN{
   #Parse and output tab-separated:
   FS="\t";
   OFS=FS;
   #Make sure array iteration order is consistent:
   PROCINFO["sorted_in"]="@ind_str_asc";
}
#First file is the .fai so that we can ensure that
# scaffolds without any BUSCOs are not purged:
FNR==NR{
   #Initialize the list of BUSCO hit pairs with all scaffold pairs:
   #This is necessary for output order and completeness.
   for (s in scaffolds) {
      if ($1 <= s) {
         loc_pairs[$1,s]=0;
      } else {
         loc_pairs[s,$1]=0;
      };
   };
   #Initialize the BUSCO hit count with all scaffolds:
   scaffolds[$1]=0;
   #Also store lengths:
   scaflens[$1]=$2;
}
#For non-commented lines, scan for Duplicated BUSCOs:
FNR<NR&&!/^#/&&$2!="Missing"{
#!/^#/&&$2!="Missing"{
   #Enter the BUSCO into a list of BUSCOs per scaffold so
   # we can easily output a list of BUSCOs to deduplicate.
   #Also add it to the appropriate count for each scaffold:
   if (length(scaffold_BUSCOs[$3]) == 0) {
      scaffold_BUSCOs[$3]=$1;
   } else {
      scaffold_BUSCOs[$3]=scaffold_BUSCOs[$3]","$1;
   };
   scaffold_composition[$3,$2]+=1;
   scaffolds[$3]+=1;
   #Keep track of the list of hit locations for each BUSCO:
   if (length(busco_loc[$1]) == 0) {
      busco_loc[$1]=$3;
   } else {
      #Check the new location against previous locations,
      # and record non-self pairs:
      #Note: If a BUSCO is found >2 times, it will contribute
      # n-m-1 times to pair counts (where m is the number of
      # self-pairs, n is the number of hits). A self-pair is
      # when two hits for the same BUSCO are found on the
      # same contig/scaffold.
      n=split(busco_loc[$1], locs, ",");
      for (l in locs) {
         if (locs[l] != $3) {
            #Only record the pair in one orientation (i.e. don't double-count):
            #Using lexicographic order here
            if (locs[l] <= $3) {
               loc_pairs[locs[l],$3]+=1;
            } else {
               loc_pairs[$3,locs[l]]+=1;
            };
         };
      };
      #Add the scaffold to the list:
      busco_loc[$1]=busco_loc[$1]","$3;
   }
}
END{
   print "ScaffoldA", "DupCountA", "BUSCOCountA", "DupRateA", "ScaffoldB", "DupCountB", "BUSCOCountB", "DupRateB", "SharedDupCount", "SharedDupRateA", "SharedDupRateB", "LengthA", "LengthB";
   for (p in loc_pairs) {
      n=split(p, locs, SUBSEP);
      #We do a lot of divide-by-zero checking below, but it really shouldn't
      # happen, since totalA or totalB == 0 implies dupsA or dupsB == 0,
      # respectively, and loc_pair[p] must also be 0. I'm just being paranoid.
      dupsA=0;
      if (((locs[1],"Duplicated") in scaffold_composition)) {
         dupsA=scaffold_composition[locs[1],"Duplicated"];
      };
      totalA=0;
      if ((locs[1] in scaffolds)) {
         totalA=scaffolds[locs[1]];
      };
      dupsB=0;
      if (((locs[2],"Duplicated") in scaffold_composition)) {
         dupsB=scaffold_composition[locs[2],"Duplicated"];
      };
      totalB=0;
      if ((locs[2] in scaffolds)) {
         totalB=scaffolds[locs[2]];
      };
      if (totalA == 0) {
         duprateA="NA";
         sharedduprateA="NA";
      } else {
         duprateA=dupsA/totalA;
         sharedduprateA=loc_pairs[p]/totalA;
      };
      if (totalB == 0) {
         sharedduprateB="NA";
         duprateB="NA";
      } else {
         duprateB=dupsB/totalB;
         sharedduprateB=loc_pairs[p]/totalB;
      };
      print locs[1], dupsA, totalA, duprateA, locs[2], dupsB, totalB, duprateB, loc_pairs[p], sharedduprateA, sharedduprateB, scaflens[locs[1]], scaflens[locs[2]];
   };
   #Now print the list of BUSCOs per contig/scaffold to STDERR:
   for (i in scaffolds) {
      if (i in scaffold_BUSCOs) {
         print i, scaffold_BUSCOs[i] > "/dev/stderr";
      } else {
         print i, "" > "/dev/stderr";
      };
   };
}
