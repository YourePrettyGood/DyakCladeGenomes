#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

#########################################################################################
# breakMisjoins.pl                                                                      #
# Version 1.0 (2019/11/12)                                                              #
# Version 1.1 (2019/11/19) Checks BED for sortedness within contig                      #
#                                                                                       #
# Description:                                                                          #
# This script breaks contigs or scaffolds apart into constituents (_ctg#) and           #
# misjoin regions (_misjoin#) based on a BED of misjoin intervals.                      #
# There's nothing too special about this script.                                        #
#                                                                                       #
# Usage:                                                                                #
#  breakMisjoins.pl [-i input_contigs.fasta] [-b misjoins.bed]                          #
# Options:                                                                              #
#  --help,-h,-?          Display this help documentation                                 #
#  --input_file,-i:      Input contigs FASTA file name (default: STDIN)                  #
#  --misjoin_bed,-b      BED file of putative misjoin intervals                          #
#  --version,-v          Output version string                                          #
#########################################################################################

my $SCRIPTNAME = "breakMisjoins.pl";
my $VERSION = "1.1";

=pod

=head1 NAME

breakMisjoins.pl - Break scaffolds or contigs at misjoin intervals provided

=head1 SYNOPSIS

breakMisjoins.pl [options]

 Options:
  --help,-h,-?          Display this help documentation
  --input_file,-i       Input contigs FASTA file name (default: STDIN)
  --misjoin_bed,-b      Sorted BED file of putative misjoin intervals
  --version,-v          Output version string

=head1 DESCRIPTION
This script breaks contigs or scaffolds apart into constituents (_ctg[number]) and
misjoin regions (_misjoin[number]) based on a BED of misjoin intervals.
There's nothing too special about this script.

=cut

my $display_version = 0;
my $help = 0;
my $man = 0;
my $input_path = "STDIN";
my $bed_path = "";
my $dispversion = 0;
my $debug = 0;
GetOptions('input_file|i=s' => \$input_path, 'misjoin_bed|b=s' => \$bed_path, 'debug|d+' => \$debug, 'help|h|?+' => \$help, man => \$man, 'version|v' => \$dispversion) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

my $contigsfh;
if ($input_path ne "STDIN") {
   unless(open($contigsfh, "<", $input_path)) {
      print STDERR "Error opening input contigs FASTA file ${input_path}.\n";
      exit 2;
   }
} else {
   open($contigsfh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $contigsfh so we can seamlessly handle piping
}

my $bedfh;
unless(open($bedfh, "<", $bed_path)) {
   print STDERR "Error opening misjoin BED file ${bed_path}.\n";
   exit 3;
}

sub breakContig {
   my $header = shift @_;
   my $sequence = shift @_;
   my @misjoins = @{shift @_};

   my $output_str = "";
   my $seq_end = length($sequence)-1;

   my $num_intervals = scalar(@misjoins);
   $output_str = ">${header}\n${sequence}\n" unless $num_intervals > 0;

   my @header_parts = split /\s+/, $header;
   my $header_prefix = $header_parts[0];

   my $good_start = 0; #Keep track of the start of the preceding non-misjoin interval
   my $i = 0;
   while ($i < $num_intervals) {
      #Only output preceding non-misjoin interval if it exists:
      if ($good_start < $misjoins[$i][0]) {
         $header_parts[0] = ">${header_prefix}_ctg" . ($i+1);
         $output_str .= join(" ", @header_parts) . "\n";
         $output_str .= substr($sequence, $good_start, $misjoins[$i][0]-$good_start) . "\n";
      }
      #Output the misjoin interval:
      $header_parts[0] = ">${header_prefix}_misjoin" . ($i+1);
      $output_str .= join(" ", @header_parts) . "\n";
      $output_str .= substr($sequence, $misjoins[$i][0], $misjoins[$i][1]) . "\n";
      #Move the good_start pointer forward to the next start:
      $good_start = $misjoins[$i][0] + $misjoins[$i][1];
      $i++;
   }
   #Catch a trailing non-misjoin interval if it exists:
   if ($good_start < $seq_end) {
      $header_parts[0] = ">${header_prefix}_ctg" . ($i+1);
      $output_str .= join(" ", @header_parts) . "\n";
      $output_str .= substr($sequence, $good_start, $seq_end - $good_start + 1) . "\n";
   }

   return $output_str;
}

my %misjoins = (); #Hash for holding misjoin intervals
print STDERR "Reading in misjoin BED ${bed_path}.\n" if $debug;
while (my $bedline = <$bedfh>) {
   chomp $bedline;
   my @bed_elems = split /\t/, $bedline;
   my ($scaf, $bedstart, $bedend) = @bed_elems[0 .. 2];
   #Check BED for sortedness within a scaffold by just asserting
   # increasing order as it's read in:
   #Basically, if you just run sort -k1,1 -k2,2n -k3,3n on your input BED,
   # you'll be fine.
   if (exists($misjoins{$scaf})) {
      my @prev_intervals = @{$misjoins{$scaf}};
      #If starts are in decreasing order, or if starts match and ends are in decreasing order:
      if ($bedstart < $prev_intervals[$#prev_intervals][0] or ($bedstart == $prev_intervals[$#prev_intervals][0] and $bedend-$bedstart < $prev_intervals[$#prev_intervals][1])) {
         print STDERR "BED intervals found out of order for ${scaf}, exiting\n";
         close($bedfh);
         exit 5;
      }
   }
   $misjoins{$scaf} = [] unless exists($misjoins{$scaf});
   push @{$misjoins{$scaf}}, [$bedstart, $bedend-$bedstart];
}
close($bedfh);

my ($header, $sequence) = ('', '');
while (my $line = <$contigsfh>) {
   chomp $line;
   #Put the contig in the hash if we've reached the next contig:
   if ($header ne '' and $line =~ />/) {
      #Only treat the part before the space as the contig/scaffold name in the BED:
      my @header_parts = split /\s+/, $header;
      print STDOUT breakContig($header, $sequence, $misjoins{$header_parts[0]}) if exists($misjoins{$header_parts[0]});
      #Feed through if no misjoins on that contig/scaffold:
      print STDOUT ">", $header, "\n", $sequence, "\n" unless exists($misjoins{$header_parts[0]});
      ($header, $sequence) = ('', '');
   }
   #Fill up the header and sequence variables if we haven't reached the next contig:
   if ($line =~ />/) {
      $header = substr $line, 1;
   } else {
      $sequence .= $line;
   }
}
#Add the final contig into the hash, if it exists:
if ($header ne '' and $sequence ne '') {
   #Only treat the part before the space as the contig/scaffold name in the BED:
   my @header_parts = split /\s+/, $header;
   print STDOUT breakContig($header, $sequence, $misjoins{$header_parts[0]}) if exists($misjoins{$header_parts[0]});
   #Feed through if no misjoins on that contig/scaffold:
   print STDOUT ">", $header, "\n", $sequence, "\n" unless exists($misjoins{$header_parts[0]});
   ($header, $sequence) = ('', '');
}
#Close the input file if it was indeed opened:
if ($input_path ne "STDIN") {
   close($contigsfh);
}

exit 0;
