#!/bin/awk -f
#Assumes unwrapped FASTA as input
#Read in the header lines:
/^>/{
   header=$0;
}
#Read in the sequence lines, store them in a hash:
!/^>/{
   seq=$0;
   seqs[header]=seq;
}
#Sort by header, and output records accordingly:
END{
   n=asorti(seqs, headers);
   for (i=1; i<=n; i++) {
      print headers[i];
      print seqs[headers[i]];
   };
}
