#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

##################################################################
#                                                                #
# Version 1.0 (2019/04/12) Initial version                       #
# Version 1.1 (2019/04/17) Assign homology based on flanking     #
#                          exon length                           #
# Version 1.2 (2019/04/22) Split arbitrary FASTAs, not just      #
#                          introns, still using SCO map          #
#                          Homology determined with other script #
# Version 1.3 (2019/05/17) Drop duplicated records, don't parse  #
#                          it based on the last record it's in.  #
##################################################################

#First pass script to split FASTAs of introns based on single copy
# orthogroups. The first line of the input TSV should be a header
# indicating the species prefixes for each column's transcript IDs
# and the first column indicates the orthogroup ID (e.g. from
# OrthoFinder).
#Version 1.2 is adapted to split any FASTAs based on header matching
# to IDs from a map TSV file. This file is generated by some other
# script (in the case of introns, findSingleCopyOrthologIntrons.pl).

my $SCRIPTNAME = "parseFASTARecords.pl";
my $VERSION = "1.3";

=pod

=head1 NAME

parseFASTARecords.pl - Split FASTAs based on a map file of sequence IDs

=head1 SYNOPSIS

parseFASTARecords.pl [options]

 Options:
  --help,-h,-?           Display this help documentation
  --fasta_list,-f        File of file names of input FASTAs (default: STDIN)
  --map,-m               Path to map TSV (e.g. of single-copy orthologs)
  --version,-v           Output version string
  --debug,-d             Output debugging info to STDERR

=head1 DESCRIPTION

This script parses FASTA records into their appropriate group-specific
FASTA based on the mapping/grouping provided. The mapping is expected
to consist of a group (the first column, e.g. an orthogroup, or an
orthologous intron ID), and record IDs that belong in that group. The
headers for the columns should be "OG" for the first column, and the
prefix used for each column (e.g. a species ID) which is separated by
an underscore from the record ID.
That sounds kind of confusing, so perhaps an example will help. Say
we have some orthogroup OG0001462 which consists of the transcripts
FBtr0340054 in species Dmel, jg6131.t1 in species Dsan_STOCAGO1482,
jg2474.t1 in species Dtei_GT53w, and jg124.t1 in species Dyak_NY73PB.
This would be represented in a mapping file as:

OG	Dmel	Dsan_STOCAGO1482	Dtei_GT53w	Dyak_NY73PB
OG0001462	FBtr0340054	jg6131.t1	jg2474.t1	jg124.t1

The input FASTAs would be expected to have headers like so:
Dmel FASTA:
>Dmel_FBtr0340054 [anything else after a whitespace is ignored]
ACCTATATGGAATATGATAGA

Dsan FASTA:
>Dsan_STOCAGO1482_jg6131.t1
ACCTACCTGGACAATGGTAGA

Dtei FASTA:
>Dtei_GT53w_jg2474.t1
ACCTACCTGTACAAGGGTAGA

Dyak FASTA:
>Dyak_NY73PB_jg124.t1
ACCTACCTGGAGAATGGTAGA

The FASTAs are provided as a file containing a newline-separated
list of paths to FASTAs.

This script opens and closes a lot of file handles, so has some
overhead there, but avoids the overhead of searching each record
of a FASTA across a list of search patterns, so should be faster
to partition FASTAs by these groups.

=cut

my $help = 0;
my $man = 0;
my $fasta_list_path = "STDIN";
my $map_path = "";
my $dispversion = 0;
my $debug = 0;
GetOptions('fasta_list|f=s' => \$fasta_list_path, 'map|m=s' => \$map_path, 'version|v' => \$dispversion, 'debug|d+' => \$debug, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

my $fasta_list_fh;
#Open the list of FASTA files, or set it up to be read from STDIN:
if ($fasta_list_path ne "STDIN") {
   unless(open($fasta_list_fh, "<", $fasta_list_path)) {
      print STDERR "Error opening list of FASTA files ${fasta_list_path}.\n";
      exit 2;
   }
} else {
   open($fasta_list_fh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $fasta_list_fh so we can seamlessly handle piping
}

#Read in the FASTA list:
my @fasta_list = ();
while (my $line = <$fasta_list_fh>) {
   chomp $line;
   push @fasta_list, $line if -e $line;
   print STDERR "FASTA file ${line} in list ${fasta_list_path} does not exist.\n" unless -e $line;
}
close($fasta_list_fh);

#Open the group map file:
my $map_fh;
unless(open($map_fh, "<", $map_path)) {
   print STDERR "Error opening map TSV ${map_path}.\n";
   exit 3;
}

#Read in the group map:
my %group_map = (); #Hash of hashes ({spp}{record id})
my %skip_IDs = (); #Hash of record IDs to omit from parsing -- meant for duplicates
my $header_line = <$map_fh>;
chomp $header_line;
my @strains = split /\t/, $header_line; #1-based indexing, as index 0 is OG
shift @strains;
my @species = map { my @elems = split /_/, $_; $elems[0] } @strains;
my $num_species = scalar(@species);
while (my $line = <$map_fh>) {
   chomp $line;
   my @records = split /\t/, $line;
   for (my $i = 1; $i <= $num_species; $i++) {
      if ($records[$i] eq "") {
         print STDERR "Species ", $species[$i-1], " is missing a record for group ", $records[0], ", cannot proceed\n";
         exit 4;
      }
      my $spp = $species[$i-1];
      my $record = $records[$i];
      my $group = $records[0];
      $group_map{$spp} = {} unless exists($group_map{$spp});
      if (exists($group_map{$spp}{$record})) {
         print STDERR "Record ${record} found multiple times for species ${spp} in group map ${map_path}.\n";
         $skip_IDs{$spp} = {} unless exists($skip_IDs{$spp});
         $skip_IDs{$spp}{$record} = undef;
         $group_map{$spp}{$record} = undef;
      }
      print STDERR "Record ${record} for species ${spp} => was ", exists($group_map{$spp}{$record}) && defined($group_map{$spp}{$record}) ? $group_map{$spp}{$record} : "", " is now ", $records[0], "\n" if $debug > 2;
      $group_map{$spp}{$record} = $group unless exists($group_map{$spp}{$record}) and defined($group_map{$spp}{$record}); #Map the record to the group
   }
}
close($map_fh);

#Drop the records that are duplicated:
for my $spp (keys %skip_IDs) {
   for my $record (keys %{$skip_IDs{$spp}}) {
      delete $group_map{$spp}{$record};
   }
}

#Now iterate through the FASTA files and sort each record into
# its appropriate file based on the map file:
#We choose to open and close file handles for each FASTA record
# since this obviates any concerns about open file handle limits.
#It may, however, lead to slower performance due to fh overhead.
for my $fasta_path (@fasta_list) {
   my $in_fasta;
   unless (open($in_fasta, "<", $fasta_path)) {
      print STDERR "Unable to open input FASTA ${fasta_path}.\n";
      exit 5;
   }
   my $fasta_header = "";
   my $fasta_seq = "";
   while (my $line = <$in_fasta>) {
      chomp $line;
      if ($line =~ /^>/) { #Header line
         if ($fasta_header ne "") {
            #Process the prior record:
            my @header_elems = split /\s+/, substr($fasta_header, 1);
            my $record_ID = $header_elems[0];
            my @recordid_elems = split /_/, $record_ID;
            print STDERR "Unable to find species for record ${record_ID}\n" unless scalar(@recordid_elems) > 1;
            my $short_recordid = pop @recordid_elems;
            #my $species = join("_", @recordid_elems);
            my $species = shift @recordid_elems;
            unless (exists($group_map{$species}{$short_recordid}) and defined($group_map{$species}{$short_recordid})) {
               print STDERR "No group found for ${record_ID}, skipping\n" if $debug > 1;
               $fasta_seq = "";
               $fasta_header = $line;
               next;
            }
            my $group = $group_map{$species}{$short_recordid};
            my $output_fn = "${group}_unwrapped.fasta";
            print STDERR "${record_ID} goes in ${output_fn}\n" if $debug > 2;
            my $output_fh;
            unless(open($output_fh, ">>", $output_fn)) {
               print STDERR "Unable to open output file ${output_fn} for appending\n";
               print STDERR "Skipping output of ${record_ID} due to append open error\n" if $debug > 2;
            } else {
               print $output_fh $fasta_header, "\n";
               print $output_fh $fasta_seq, "\n";
               close($output_fh);
            }
         }
         $fasta_seq = "";
         $fasta_header = $line;
      } else { #Sequence line
         $fasta_seq .= $line;
      }
   }
   #Process the last record:
   my @header_elems = split /\s+/, substr($fasta_header, 1);
   my $record_ID = $header_elems[0];
   my @recordid_elems = split /_/, $record_ID;
   my $short_recordid = pop @recordid_elems;
   #my $species = join("_", @recordid_elems);
   my $species = shift @recordid_elems;
   if (exists($group_map{$species}{$short_recordid}) and defined($group_map{$species}{$short_recordid})) {
      my $group = $group_map{$species}{$short_recordid};
      my $output_fn = "${group}_unwrapped.fasta";
      print STDERR "${record_ID} goes in ${output_fn}\n" if $debug > 2;
      my $output_fh;
      unless(open($output_fh, ">>", $output_fn)) {
         print STDERR "Unable to open output file ${output_fn} for appending\n";
         print STDERR "Skipping output of ${record_ID} due to append open error\n" if $debug > 2;
      } else {
         print $output_fh $fasta_header, "\n";
         print $output_fh $fasta_seq, "\n";
         close($output_fh);
      }
   } else {
      print STDERR "No group found for ${record_ID}, skipping\n" if $debug > 1;
   }
   close($in_fasta);
}

exit 0;
