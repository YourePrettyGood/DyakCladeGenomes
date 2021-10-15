#!/bin/awk -f
#spid is the only input variable, and it's the species ID
# (the ID of the reference genome, so e.g. Dtei_GT53w for Dtei)
BEGIN{
   FS="\t";
   OFS=FS;
}
{
   print $1, $2-1, $3, "OGID="$5";"spid"ID="$6":"$7, ".", $4;
}
