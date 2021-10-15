#!/bin/awk -f
#Generate simple summaries of a FASTQ file, like the total number of reads,
# average read length, and total base count
BEGIN{
   OFS="\t";
}
NR%4==2{
   sum+=length($0);
   count+=1;
}
END{
   print count, sum/count, sum;
}
