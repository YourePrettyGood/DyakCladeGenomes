#!/bin/awk -f
#Unwraps a protein FASTA and adds a stop codon (*) to the end
# of each protein sequence.
#num is used to keep track of original sequence order
BEGIN{
   num=1;
}
#Use a simplified header (first space-separated token without >)
# as the array key, and store the header and it's order:
/^>/{
   header=substr($1, 2);
   seqorder[num]=header;
   headers[header]=$0;
   num+=1;
}
#Store the protein sequence progressively, implicitly unwrapping:
!/^>/{
   seqs[header]=seqs[header]""$0;
}
#Output the protein sequences in the input order:
END{
   for (i=1; i < num; i++) {
      print headers[seqorder[i]];
      print seqs[seqorder[i]]"*";
   };
}
