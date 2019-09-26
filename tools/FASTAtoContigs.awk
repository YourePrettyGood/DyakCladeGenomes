#!/bin/awk -f
/^>/{
   header=$0;
}
!/^>/{
#For simplicity here, all we do is split on a string of Ns of length >= 1
# although we don't account for lowercase n
#Fancier versions could use {#} instead of + for splitting on fixed-length
# gaps, or {#,#} for variable-length gaps, and more complex regexes for
# gaps with varying case and/or character sets.
#For instance, to split only on what NCBI considers "unknown" gaps (100 Ns):
# split($0, ctgs, /[N]{100}/);
   split($0, ctgs, /[N]+/);
   for (ctg in ctgs) {
      print header"_ctg"ctg;
      print ctgs[ctg];
   };
}
