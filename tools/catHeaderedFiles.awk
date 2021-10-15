#!/bin/awk -f
#This script does a very simple but oft-overlooked task of concatenating
# headered files. Basically, the whole first file is printed out, and
# then all but the first line of every subsequnt file is printed.
#I don't know of any base Unix command that does this -- please let me
# know if there's an obvious one that works for > 2 files!
FNR==NR{
   print;
}
FNR<NR&&FNR>1{
   print;
}
