#!/bin/awk -f
/^>/{
   skip=0;
   split(substr($0, 2), idparts, "_");
#Slightly complicated regex, but basically just covers CY*, NY*, and R9b
   if (idparts[1] == "Dyak" && idparts[2] ~ /^[RCN][Y9][0-9b]+[A-C]?[0-9]?/) {
      skip=1;
#And an extra regex for the other inbred Dyak lines (Jess, ST, 2.8.1)
   } else if (idparts[1] == "Dyak" && idparts[2] ~ /^(Jess|ST|2\.8\.1)/) {
      skip=1;
#Also covers Dsim MD* and NS* (2019/06/25 edit)
#Might be slightly overzealous for generality, but for our data it's fine
   } else if (idparts[1] == "Dsim" && idparts[2] ~ /^[MN][DS][0-9]+/) {
      skip=1;
#Added Dmel DPGP3 haploid embryo sequences (2019/11/17 edit)
# (Not a full regex match, but at least covers ZI)
   } else if (idparts[1] == "Dmel" && idparts[2] ~ /^ZI/) {
      skip=1;
#Added Dtei TUZ11 (2021/05/30 edit)
   } else if (idparts[1] == "Dtei" && idparts[2] == "TUZ11") {
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
