#!/bin/awk -f
BEGIN{
   gsub(/_/, " ", qid);
   gsub(/_/, " ", rid);
}
NR==1{
   print "set terminal png size 1400,1400 font \"arial,24\"";
}
NR>1{
   if ($0 ~ "QRY") {
      print "set ylabel \""qid"\"";
   } else if ($0 ~ "REF") {
      print "set xlabel \""rid"\"";
   } else if ($0 ~ "set [xy]tics ") {
      fixscafnames=1;
      print;
   } else if ($0 ~ ")") {
      fixscafnames=0;
      print;
   } else {
      if (fixscafnames>0) {
         gsub(/_/, " ", $0);
      }
      print $0;
   }
}
