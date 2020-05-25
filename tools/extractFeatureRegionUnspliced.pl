#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
# Version 1.0 (2019/12/08) Derived from extractRegions.pl      #
################################################################

#First pass script to extract features of a given type in a GFF3 file
# from a FASTA file as individual records with headers based on the ID.

my $SCRIPTNAME = "extractFeatureRegionsUnspliced.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

extractFeatureRegionsUnspliced.pl - Extract GFF3 feature regions from a FASTA

=head1 SYNOPSIS

extractFeatureRegionsUnspliced.pl [options]

 Options:
  --help,-h,-?          Print this help documentation
  --input_genome,-i     Path to input genome FASTA file (default: STDIN)
  --gff3_path,-g        Path to GFF3 file indicating regions to be extracted
                        from the FASTA
  --feature,-f          Feature type (column 3 of GFF3) to extract
  --prefix,-p           Prefix for FASTA headers (off by default)
  --version,-v          Output version string

=head1 DESCRIPTION

This script extracts the features at positions given by the input GFF3 file
from the input FASTA file, and outputs each feature as its own record in a
FASTA format output to STDOUT.  Note that splicing and parental relationships
are ignored, so extracting an mRNA record will include introns.

=cut

sub revcomp($) {
   my $input_sequence = shift @_;
   my $reverse_sequence = reverse $input_sequence; #Reverse
   $reverse_sequence =~ tr/AaCcGgTtRrYySsWwKkMmBbDdHhVvNn/TtGgCcAaYyRrSsWwMmKkVvHhDdBbNn/; #Complement incl. IUPAC degenerate bases
   return $reverse_sequence;
}

my $help = 0;
my $man = 0;
my $genome_path = "STDIN";
my $gff3_path = "";
my $feature_type = "";
my $header_prefix = "";
my $dispversion = 0;
GetOptions('input_genome|i=s' => \$genome_path, 'gff3_file|g=s' => \$gff3_path, 'feature|f=s' => \$feature_type, 'prefix|p=s' => \$header_prefix, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;
$header_prefix .= "_" unless $header_prefix eq "";

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

#Open the GFF3 file:
my $gff3_fh;
unless(open($gff3_fh, "<", $gff3_path)) {
   print STDERR "Error opening GFF3 file ${gff3_path}.\n";
   exit 2;
}

my $FASTA_skip = 0; #Skip all lines after the ##FASTA line if present
my $feature_index = 1; #Use this feature index to construct an ID if ID unavailable
my $prev_Parent = ""; #Use this to reset feature_index if we change Parents
my %regions_per_scaffold = ();
while (my $line = <$gff3_fh>) {
   chomp $line;
   $FASTA_skip = 1 if $line =~ /^##FASTA/; #Trigger FASTA skip
   next if $FASTA_skip;
   next if $line =~ /^#/; #Skip header lines
   my ($scaffold, $source, $feature, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $line, 9;
   #Skip any features that aren't of the desired type:
   next unless $feature eq $feature_type;
   #Grab the feature ID if available, or the parent if ID isn't available:
   my ($ID, $Parent);
   my @tags = split /;[ ]*/, $attributes;
   for my $tag (@tags) {
      my ($tag_name, $tag_value) = split /=/, $tag, 2;
      if ($tag_name eq "ID") {
         $ID = $tag_value;
      } elsif ($tag_name eq "Parent") {
         $Parent = $tag_value;
      }
   }
   #Infer ID if ID isn't in the attributes string:
   unless (defined($ID)) {
      $feature_index = 1 unless $Parent eq $prev_Parent;
      $ID = "${Parent}.${feature_type}${feature_index}";
      $feature_index++;
      $prev_Parent = $Parent;
   }
   my ($start_0_based, $end_1_based) = ($start-1, $end);
   push @{$regions_per_scaffold{$scaffold}}, "${ID}=${start_0_based}:".($end_1_based-$start_0_based)."=${strand}";
}

close($gff3_fh);

#Now we can iterate through the genome FASTA, and output regions as they arise
my $scaffold_name = "";
my $scaffold_sequence = "";
while (my $line = <$genome_fh>) {
   chomp $line;
   #If we're at a header line and we've seen header lines before,
   # output the sites from the previous scaffold (since we're on
   # a new scaffold's header line):
   if ($line =~ /^>/) {
      if ($scaffold_name ne "" and exists($regions_per_scaffold{$scaffold_name})) {
         #Output the regions on this scaffold:
         for my $locus (@{$regions_per_scaffold{$scaffold_name}}) {
            my ($locus_name, $range, $strand) = split /=/, $locus, 3;
            my ($range_start, $range_length) = split /:/, $range, 2;
            my $region_sequence = substr $scaffold_sequence, $range_start, $range_length;
            $region_sequence = revcomp($region_sequence) if $strand eq "-";
            my $region_name = "";
            $region_name = ">${header_prefix}";
            $region_name .= "${locus_name} ${scaffold_name}:" . ($range_start+1) . ".." . ($range_start+$range_length) . "(${strand})";
            print STDOUT $region_name, "\n", $region_sequence, "\n";
         }
      }
      my $scaffold_name_line = substr $line, 1; #Get rid of the prefixed ">"
      #Only preseve the first word of the scaffold header:
      my @scaffold_words = split /\s+/, $scaffold_name_line;
      $scaffold_name = $scaffold_words[0];
      $scaffold_sequence = ""; #Clear out the old sequence
   } else { #Sequence line
      $scaffold_sequence .= $line;
   }
}
close($genome_fh);

#Now make sure we account for the last scaffold:
if (exists($regions_per_scaffold{$scaffold_name})) {
   #Concatenate the sites:
   for my $locus (@{$regions_per_scaffold{$scaffold_name}}) {
      my ($locus_name, $range, $strand) = split /=/, $locus, 3;
      my ($range_start, $range_length) = split /:/, $range, 2;
      my $region_sequence = substr $scaffold_sequence, $range_start, $range_length;
      $region_sequence = revcomp($region_sequence) if $strand eq "-";
      my $region_name = "";
      $region_name = ">${header_prefix}";
      $region_name .= "${locus_name} ${scaffold_name}:" . ($range_start+1) . ".." . ($range_start+$range_length) . "(${strand})";
      print STDOUT $region_name, "\n", $region_sequence, "\n";
   }
}

exit 0;
