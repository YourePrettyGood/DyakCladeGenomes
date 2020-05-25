#!/bin/awk -f
BEGIN{
   FS="\t";
   OFS=FS;
}
#First file is the reference GFF3 containing chromosome_band records:
FNR==NR&&$3=="chromosome_band"{
   split($9, tags, ";");
   for (i in tags) {
      split(tags[i], elems, "=");
      if (elems[1] == "Name") {
         split(elems[2], bandname, "-");
         bands[$1":"($4-1)"-"$5]=bandname[2];
#         print $1":"($4-1)"-"$5"\t"bandname[2] > "/dev/stderr";
      };
   };
}
#Second file is the BED6 of SCOs for the reference:
FNR<NR{
   bandstr="";
   delete band_hits;
   hitcount=0;
   #Add a band to the list if any part of the SCO overlaps with the band:
   for (interval in bands) {
      split(interval, scafpart, ":");
      if ($1 == scafpart[1]) {
         split(scafpart[2], intparts, "-");
         if (($2 >= intparts[1] && $2 < intparts[2]) || ($3 <= intparts[2] && $3 > intparts[1])) {
            band_hits[++hitcount]=bands[interval];
#            print $1"\t"scafpart[1]"\t"$2"\t"intparts[1]"\t"$3"\t"intparts[2]"\t"hitcount"\t"interval"\t"bands[interval] > "/dev/stderr";
         };
      };
   };
   #Make sure to sort the bands before outputting:
   n_hits=asort(band_hits);
   for (i=1; i<=n_hits; i++) {
      bandstr=bandstr","band_hits[i];
   };
   #Append list of bands to the Name column:
   print $1, $2, $3, $4";bands="substr(bandstr, 2), $5, $6;
}
