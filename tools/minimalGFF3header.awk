#!/bin/awk -f
#A minimal GFF3 header includes the gff version line, and
# sequence-region headers
BEGIN{
   print "##gff-version 3";
}
{
   print "##sequence-region "$1" 1 "$2;
}
