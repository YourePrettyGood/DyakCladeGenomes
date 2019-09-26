#!/bin/awk -f
/^>/{
   skip=0;
   split(substr($0, 2), idparts, "_");
#Slightly complicated regex, but basically just covers CY*, NY*, and R9b
#Also covers Dsim MD* and NS* (2019/06/25 edit)
#Might be slightly overzealous for generality, but for our data it's fine
   if (idparts[1] == "Dyak" && idparts[2] ~ /^[RCN][Y9][0-9b]+[A-C]?[0-9]?/) {
      skip=1;
   } else if (idparts[1] == "Dsim" && idparts[2] ~ /^[MN][DS][0-9]+/) {
      skip=1;
   } else {
      print;
   };
}
!/^>/{
   if (skip == 0) {
      print;
   };
}
