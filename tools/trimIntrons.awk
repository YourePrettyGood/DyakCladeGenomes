#!/bin/awk -f
BEGIN{
#Default to trimming 2 bp from each side:
   if (length(ltrim) == 0) {
      ltrim=2;
   };
   if (length(rtrim) == 0) {
      rtrim=2;
   };
}
/^>/{
#Could modify this so that we append something to the header to
# indicate that it has been trimmed:
   print;
}
!/^>/{
#Trim ltrim off left end, rtrim off right end,
# and convert all bases to uppercase:
   trimmed=toupper(substr($0, 1+ltrim, length($0)-ltrim-rtrim));
   print trimmed;
}
