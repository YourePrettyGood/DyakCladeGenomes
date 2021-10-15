#!/bin/awk -f
#This script is a bit customized for Drosophila arm nomenclature,
# but the idea is to generate a table of offsets for major
# chromosome arms between species based on the syntenic single-copy
# orthogroup closest to the telomere (so smallest coordinate for
# X, 2L, 3L, and 4, but the largest coordinate for 2R and 3R).
#We also do a bit of customized trimming of strings, like getting
# rid of the pesky Scf_ prefix for all scaffold names in the Dsim
# FlyBase assembly, and removing the _proteome suffix from the
# assembly ID column of orthogroup_locations.
#We determine left- versus right-alignment of scaffolds based on
# the presence of an R at the end of the scaffold ID.
#
#The input file must be scaffold- and position-sorted, and must only
# contain single-copy orthogroups. This script does not check either
# characteristic, though the Makefile in the genomewide_pi directory
# should produce a file with only SCOs (i.e. orthogroup_locations),
# and the sorting requirement only involves a command like this:
# sort -k1,1 -k2,2n -k3,3n < [path to]/orthogroup_locations
BEGIN{
   FS="\t";
   OFS=FS;
   nspp=0;
}
#Columns of the orthogroup_locations file:
# 1) Scaffold ID
# 2) Leftmost CDS position in GFF3 (1-based)
# 3) Rightmost CDS position in GFF3 (1-based)
# 4) Strand of the CDS in GFF3
# 5) Orthogroup ID (i.e. from OrthoFinder)
# 6) Assembly ID (typically [spp]_[strain]_proteome)
# 7) Transcript ID as used by OrthoFinder
{
   #Start by trimming what needs to be trimmed:
   sub(/^Scf_/, "", $1);
   sub(/_proteome$/, "", $6);
   #Now establish the map of the scaffolds and their handedness:
   if ($1 ~ /R$/) {
      arms[$1]="R";
   } else {
      arms[$1]="L";
   };
   #Fill in the orthogroup-to-scaffold map:
   if ($5 in orthogroups) {
      #This covers the case where the SCO isn't syntenic, so we simply
      # append the other scaffold to make a list. The only point of the
      # list is to serve as a filter later for syntenic SCOs:
      if ($1 != orthogroups[$5]) {
         orthogroups[$5]=orthogroups[$5]","$1;
      };
   } else {
      orthogroups[$5]=$1;
   };
   #Keep track of the species/assemblies in the input file for later iteration:
   if (!($6 in spp)) {
      spp[$6]=++nspp;
   };
   #Store the arm- and species/assembly-specific rank of this orthogroup
   # (where the rank is determined based on the left of the scaffold):
   #We also simultaneously keep track of the maximum per-arm per-assembly
   # rank, as this will be used for calculating the right-based rank.
   pos[$6,$5]=++maxrank[$6,$1];
   #Store the left-most position of the CDS for that assembly and orthogroup:
   #Using left-most versus right-most position shouldn't matter much.
   loc[$6,$5]=$2;
}
END{
   PROCINFO["sorted_in"]="@ind_str_asc";
   #Iterate over the orthogroups to only retain syntenic orthogroups and
   # then calculate their rank sums, finding the orthogroup on each arm
   # with the minimum rank sum (i.e. the left- or right-most syntenic SCO):
   for (og in orthogroups) {
      #Immediately skip any non-syntenic SCOs:
      syntenic=split(orthogroups[og], ogscafs, ",");
      if (syntenic != 1) {
         continue;
      };
      #Calculate the sum of ranks for the given orthogroup across assemblies:
      ranksum=0;
      if (arms[ogscafs[1]] == "R") {
         for (sp in spp) {
            #We invert the left-most rank using the max rank to get the
            # right-most rank:
            ranksum+=maxrank[sp,ogscafs[1]]-pos[sp,og]+1;
         };
      } else {
         for (sp in spp) {
            ranksum+=pos[sp,og];
         };
      };
      #Now determine the SCO on each arm that minimizes this rank sum:
      #Only update the lineup and lineupranksum hashes if no value exists,
      # or if the current rank sum is less than the stored value.
      if (!(ogscafs[1] in lineup) || lineupranksum[ogscafs[1]] > ranksum) {
         lineup[ogscafs[1]]=og;
         lineupranksum[ogscafs[1]]=ranksum;
      };
   };
   #Now determine the offsets for each assembly-arm combination by
   # finding the per-assembly per-arm maximum leftmost coordinate
   print "Species", "Assembly", "Scaffold", "Offset";
   for (arm in lineup) {
      maxstart="";
      for (sp in spp) {
         if (length(maxstart) == 0 || loc[sp,lineup[arm]] > maxstart) {
            maxstart=loc[sp,lineup[arm]];
         };
      };
      #Print out the offsets:
      for (sp in spp) {
         split(sp, spstrain, "_");
         if (length(debug) > 0) {
            print spstrain[1], sp, arm, maxstart-loc[sp,lineup[arm]], loc[sp,lineup[arm]], lineup[arm], pos[sp,lineup[arm]];
         } else {
            print spstrain[1], sp, arm, maxstart-loc[sp,lineup[arm]];
         };
      };
   };
}
