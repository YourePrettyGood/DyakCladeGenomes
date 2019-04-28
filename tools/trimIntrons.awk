#!/bin/awk -f
BEGIN{
#For the time being, we hard-code the trim_window value
#TODO: Figure out how to detect if -v passes it in, and don't overwrite
   trim_window=2;
}
/^>/{
#Could modify this so that we append something to the header to
# indicate that it has been trimmed:
   print;
}
!/^>/{
#This line simply trims one trim_window length off each end,
# and converts all bases to uppercase:
   trimmed=toupper(substr($0, 1+trim_window, length($0)-2*trim_window));
   print trimmed;
}
