#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

####################################################################
#                                                                  #
# Version 1.0 (2019/10/20) Initial version                         #
####################################################################

#Quick script to identify intervals soft-masked in an assembly
# and output these intervals as a BED.
#Order of the BED is the same as the input assembly FASTA

my $SCRIPTNAME = "identifySoftmaskedRegions.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

identifySoftmaskedRegions.pl - Generate a BED of softmasked regions in a haploid assembly FASTA

=head1 SYNOPSIS

identifySoftmaskedRegions.pl [options]

 Options:
  --help,-h,-?           Display this help documentation
  --input_genome,-i      Path to input genome FASTA file (default: STDIN)
  --version,-v           Output version string
  --debug,-d             Output debugging info to STDERR

=head1 DESCRIPTION

This script generates a BED of softmasked intervals from an input haploid genome
assembly FASTA.

=cut

my $help = 0;
my $man = 0;
my $genome_path = "STDIN";
my $dispversion = 0;
my $debug = 0;
GetOptions('input_genome|i=s' => \$genome_path, 'version|v' => \$dispversion, 'debug|d+' => \$debug, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

#Open the genome FASTA file, or set it up to be read from STDIN:
my $genome_fh;
if ($genome_path ne "STDIN") {
   unless(open($genome_fh, "<", $genome_path)) {
      print STDERR "Error opening genome FASTA file ${genome_path}.\n";
      exit 1;
   }
} else {
   open($genome_fh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $genome_fh so we can seamlessly handle piping
}

#To mitigate issues with line wrapping, we read each scaffold in, and then
# identify the softmasked regions for that scaffold. The code is slightly
# messier, but more robust to input variation.
my $scaffold_name = "";
my $scaffold_sequence = "";
while (my $line = <$genome_fh>) {
   chomp $line;
   #If we're at a header line and we've seen header lines before,
   # output the softmasked regions from the previous scaffold (since
   # we're on a new scaffold's header line):
   if ($line =~ /^>/) {
      if ($scaffold_name ne "") {
         #Search the scaffold for softmasked regions and output
         # them progressively:
         #Note: @- array is 0-based index of start of match
         # @+ array is 0-based index of one past end of match
         #So this naturally produces a BED interval
         print join("\t", $scaffold_name, $-[0], $+[0]), "\n" while $scaffold_sequence =~ /[^ACGTN]+/g;
      }
      $scaffold_name = substr $line, 1; #Get rid of the prefixed ">"
      $scaffold_sequence = ""; #Clear out the old sequence
   } else { #Sequence line
      $scaffold_sequence .= $line;
   }
}
close($genome_fh);
#Now make sure we account for the last scaffold:
if ($scaffold_name ne "") {
   #Search the scaffold for softmasked regions and output
   # them progressively:
   #Note: @- array is 0-based index of start of match
   # @+ array is 0-based index of one past end of match
   #So this naturally produces a BED interval
   print join("\t", $scaffold_name, $-[0], $+[0]), "\n" while $scaffold_sequence =~ /[^ACGTN]+/g;
}

exit 0;
