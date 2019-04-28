#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

use Data::Dumper;

################################################################
#                                                              #
# Version 1.0 (2019/04/22) Assign homology based on flanking   #
#                          exon length                         #
################################################################

#First pass script to split FASTAs of introns based on single copy
# orthogroups. The first line of the input TSV should be a header
# indicating the species prefixes for each column's transcript IDs
# and the first column indicates the orthogroup ID (e.g. from
# OrthoFinder).

my $SCRIPTNAME = "findSingleCopyOrthologIntrons.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

findSingleCopyOrthologs.pl - Identify orthologous introns based on single-copy orthogroups

=head1 SYNOPSIS

findSingleCopyOrthologIntrons.pl [options]

 Options:
  --help,-h,-?           Display this help documentation
  --fasta_list,-f        File of file names of input FASTAs (default: STDIN)
  --sco_map,-m           Path to single-copy orthogroup map TSV
  --threshold,-t         Maximum relative deviance for homology matching
                         (default: 0.10, i.e. 10%)
  --intron_len,-i        Use intron length filter in addition to flanking
                         exon length filters
  --version,-v           Output version string
  --debug,-d             Output debugging info to STDERR

=head1 DESCRIPTION

This script identifies orthologous introns based on the gene mapping
between species provided in the "single-copy orthogroup map" (e.g.
from the output of OrthoFinder after converting OrthoFinder protein
IDs back to the transcript IDs from the original annotation GFF3s).
The "single-copy orthogroup map" must have a header line consisting
of a first column labeled "OG", and subsequent columns labeled with
the species prefix used in extracting introns. Orthogroup lines after
the header consist of transcript IDs without the species prefix.
The output is a file similar to the "single-copy orthogroup map",
consisting of the same header line, and rows indicating an orthologous
intron ID (i.e. [orthogroup ID]_intron[#]), and the orthologous
introns using their within-species intron IDs. This output can then
be used with parseSingleCopyOrthologs.pl to parse introns into FASTAs
for their separate orthologous groups.

=cut

my $help = 0;
my $man = 0;
my $fasta_list_path = "STDIN";
my $sco_map_path = "";
my $max_rel_diff = 0.10;
my $intron_length_filter = 0;
my $dispversion = 0;
my $debug = 0;
GetOptions('fasta_list|f=s' => \$fasta_list_path, 'sco_map|m=s' => \$sco_map_path, 'threshold|t=f' => \$max_rel_diff, 'intron_len|i' => \$intron_length_filter, 'version|v' => \$dispversion, 'debug|d+' => \$debug, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

print STDERR "Maximum relative deviance must be non-negative, not ${max_rel_diff}\n" if $max_rel_diff < 0;
exit 1 if $max_rel_diff < 0;

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

#Open the Single-copy orthogroup map file:
my $sco_map_fh;
unless(open($sco_map_fh, "<", $sco_map_path)) {
   print STDERR "Error opening single-copy orthogroup map TSV ${sco_map_path}.\n";
   exit 3;
}

#Read in the single-copy orthogroup map:
my %orthogroup_map = (); #Maps from species-specific transcript IDs to orthogroups
my %txid_map = (); #Hash of hashes ({OG}{spp}) maps from orthogroups to species-specific transcript IDs
my $header_line = <$sco_map_fh>;
chomp $header_line;
my @SCO_species = split /\t/, $header_line; #Index 0 is OG, so remove it
shift @SCO_species;
my $num_species = scalar(@SCO_species);
while (my $line = <$sco_map_fh>) {
   chomp $line;
   my @transcripts = split /\t/, $line;
   $txid_map{$transcripts[0]} = {} unless exists($txid_map{$transcripts[0]});
   for (my $i = 1; $i <= $num_species; $i++) {
      my $spp_transcript_ID = join("_", $SCO_species[$i-1], $transcripts[$i]);
      print STDERR "Transcript ${spp_transcript_ID} found multiple times in single-copy orthogroup map ${sco_map_path}.\n" if exists($orthogroup_map{$spp_transcript_ID});
      print STDERR "Transcript ${spp_transcript_ID} => ", $transcripts[0], "\n" if $debug > 2;
      $orthogroup_map{$spp_transcript_ID} = $transcripts[0]; #Map the transcript to the orthogroup
      $txid_map{$transcripts[0]}{$SCO_species[$i-1]} = $transcripts[$i]; #Map the orthogroup to the transcript
   }
}
close($sco_map_fh);

#Keep track of the species actually found in the FASTAs (may not match SCO map:
my %species_found = ();

#Do a first pass through the intron FASTAs to collect flanking
# exon information for the basic homology scans:
my %intron_metadata = (); #Hash of hashes of arrays of hashes ({OG}{spp}[intron]{il,lel,rel}=length)
for my $fasta_path (@fasta_list) {
   my $intron_fasta;
   unless (open($intron_fasta, "<", $fasta_path)) {
      print STDERR "Unable to open intron FASTA ${fasta_path}.\n";
      exit 4;
   }
   my $fasta_header = "";
   while (my $line = <$intron_fasta>) {
      chomp $line;
      if ($line =~ /^>/) { #Header line
         if ($fasta_header ne "") {
            #Process the prior record:
            my @header_elems = split /\s+/, substr($fasta_header, 1);
            if ($header_elems[0] =~ /^(\S+)\.intron(\d+)$/) {
               my $transcript_ID = $1;
               my $intron_num = $2;
               unless (exists($orthogroup_map{$transcript_ID})) {
                  print STDERR "No orthogroup found for ${transcript_ID}, skipping\n" if $debug > 1;
                  $fasta_header = $line;
                  next;
               }
               my @txid_elems = split /_/, $header_elems[0];
               my $short_txid = pop @txid_elems;
               my $species = join("_", @txid_elems);
               $species_found{$species} = undef;
               my $OG = $orthogroup_map{$transcript_ID};
               my $intron_length = $header_elems[2];
               my $left_exon_length = $header_elems[3];
               my $right_exon_length = $header_elems[4];
               $intron_metadata{$OG} = {} unless exists($intron_metadata{$OG});
               $intron_metadata{$OG}{$species} = {} unless exists($intron_metadata{$OG}{$species});
               print STDERR "Redundant intron found for ${transcript_ID} intron ${intron_num}, overwriting\n" if exists($intron_metadata{$OG}{$species}{$intron_num});
               $intron_metadata{$OG}{$species}{$intron_num} = {'il'=>$intron_length, 'lel'=>$left_exon_length, 'rel'=>$right_exon_length} unless exists($intron_metadata{$OG}{$species}{$intron_num});
            } else {
               print STDERR "Unable to parse intron FASTA header for ${fasta_path} ", $header_elems[0], ", skipping\n" if $debug;
            }
         }
         $fasta_header = $line;
      }
   }
   #Process the last record:
   my @header_elems = split /\s+/, substr($fasta_header, 1);
   if ($header_elems[0] =~ /^(\S+)\.intron(\d+)$/) {
      my $transcript_ID = $1;
      my $intron_num = $2;
      print STDERR "No orthogroup found for ${transcript_ID}, skipping\n" if (!exists($orthogroup_map{$transcript_ID}) and $debug > 1);
      if (exists($orthogroup_map{$transcript_ID})) {
         my @txid_elems = split /_/, $header_elems[0];
         my $short_txid = pop @txid_elems;
         my $species = join("_", @txid_elems);
         $species_found{$species} = undef;
         my $OG = $orthogroup_map{$transcript_ID};
         my $intron_length = $header_elems[2];
         my $left_exon_length = $header_elems[3];
         my $right_exon_length = $header_elems[4];
         $intron_metadata{$OG} = {} unless exists($intron_metadata{$OG});
         $intron_metadata{$OG}{$species} = {} unless exists($intron_metadata{$OG}{$species});
         print STDERR "Redundant intron found for ${transcript_ID} intron ${intron_num}, overwriting\n" if exists($intron_metadata{$OG}{$species}{$intron_num});
         $intron_metadata{$OG}{$species}{$intron_num} = {'il'=>$intron_length, 'lel'=>$left_exon_length, 'rel'=>$right_exon_length} unless exists($intron_metadata{$OG}{$species}{$intron_num});
      }
   } else {
      print STDERR "Unable to parse intron FASTA header for ${fasta_path} ", $header_elems[0], ", skipping\n" if $debug;
   }
   close($intron_fasta);
}

#Function for calculating relative difference:
sub reldiff {
   my $expected = shift @_;
   my $observed = shift @_;
   return((abs($observed-$expected)+0.0)/$observed); #Make sure we force float division
}

#Run the basic homology scans:
my %intron_homology = (); #Hash of hashes of hashes ({OG}{intron ID}{species}=species-specific intron ID)
my $MAX_SENTINEL = $max_rel_diff + 1.0; #Use a large number as sentinel for failed min scan
print STDERR "What are you, crazy?  Filtering at ${max_rel_diff} relative difference is gigantic!\n" unless $max_rel_diff < 20.0;

for my $OG (keys %intron_metadata) {
   my @species = keys %{$intron_metadata{$OG}};
   my $spp1 = shift @species;
   my @spp1_introns = keys %{$intron_metadata{$OG}{$spp1}};
   my $spp1_intron_count = scalar(@spp1_introns);
   #Make the hash slot even if we don't have homology across all species:
   $intron_homology{$OG} = {} unless exists($intron_homology{$OG});
   #Identify homologous introns for each intron in species 1, should be transitive:
   for (my $h=0; $h < $spp1_intron_count; $h++) {
      #Plug in the species 1 intron into the map:
      $intron_homology{$OG}{$spp1_introns[$h]} = {} unless exists($intron_homology{$OG}{$spp1_introns[$h]});
      print STDERR "Duplicate homologous intron added for ${OG} ${spp1} intron ", $spp1_introns[$h], "\n" if $debug > 1 and exists($intron_homology{$OG}{$spp1_introns[$h]}{$spp1});
      $intron_homology{$OG}{$spp1_introns[$h]}{$spp1} = $spp1_introns[$h];
      print STDERR "Placed ${OG} ${spp1} intron ", $spp1_introns[$h], " into map, mapping to ", $spp1_introns[$h], "\n" if $debug > 3;
      #Iterate over each remaining species:
      for my $spp (@species) {
         #Iterate over each intron in a focal species:
         my @spp_introns = keys %{$intron_metadata{$OG}{$spp}};
         my $spp_intron_count = scalar(@spp_introns);
         #Initialize some variables for scanning for minimum relative differences:
         my ($min_ird, $min_lerd, $min_rerd) = ($MAX_SENTINEL, $MAX_SENTINEL, $MAX_SENTINEL);
         my ($min_ird_i, $min_lerd_i, $min_rerd_i) = (-1, -1, -1); #Index of the intron
         for (my $i = 0; $i < $spp_intron_count; $i++) {
            my $intron_reldiff = reldiff($intron_metadata{$OG}{$spp1}{$spp1_introns[$h]}{'il'}, $intron_metadata{$OG}{$spp}{$spp_introns[$i]}{'il'});
            my $left_exon_reldiff = reldiff($intron_metadata{$OG}{$spp1}{$spp1_introns[$h]}{'lel'}, $intron_metadata{$OG}{$spp}{$spp_introns[$i]}{'lel'});
            my $right_exon_reldiff = reldiff($intron_metadata{$OG}{$spp1}{$spp1_introns[$h]}{'rel'}, $intron_metadata{$OG}{$spp}{$spp_introns[$i]}{'rel'});
            print STDERR "${OG} ${spp1} intron ", $spp1_introns[$h], " ${spp} intron ", $spp_introns[$i], " ${left_exon_reldiff} ${right_exon_reldiff} ${intron_reldiff}\n" if $debug > 3;
            if ($intron_reldiff < $min_ird) { #Maybe skip this filter?
               $min_ird = $intron_reldiff;
               $min_ird_i = $i;
            }
            if ($left_exon_reldiff < $min_lerd) {
               $min_lerd = $left_exon_reldiff;
               $min_lerd_i = $i;
            }
            if ($right_exon_reldiff < $min_rerd) {
               $min_rerd = $right_exon_reldiff;
               $min_rerd_i = $i;
            }
         }
         #Apply the filters:
         print STDERR "${OG} ${spp1} intron ", $spp1_introns[$h], " ${spp} intron ", $spp_introns[$min_lerd_i], " ${min_lerd_i} ${min_rerd_i} ${min_ird_i} ${min_lerd} ${min_rerd} ${min_ird}\n" if $debug > 2;
         if ($min_lerd_i == $min_rerd_i and $min_lerd_i >= 0) {
            next if ($intron_length_filter and $min_ird_i != $min_lerd_i);
            next if ($intron_length_filter and $min_ird > $max_rel_diff);
            if ($min_lerd <= $max_rel_diff and $min_rerd <= $max_rel_diff) {
               #The following two lines should never execute given that we've already added the spp1 intron:
               $intron_homology{$OG} = {} unless exists($intron_homology{$OG});
               $intron_homology{$OG}{$spp1_introns[$h]} = {} unless exists($intron_homology{$OG}{$spp1_introns[$h]});
               print STDERR "Duplicate homologous intron added for ${OG} ${spp} intron ", $spp_introns[$min_lerd_i], "\n" if $debug > 1 and exists($intron_homology{$OG}{$spp1_introns[$h]}{$spp});
               $intron_homology{$OG}{$spp1_introns[$h]}{$spp} = $spp_introns[$min_lerd_i];
               print STDERR "Placed ${OG} ${spp} intron ", $spp_introns[$min_lerd_i], " into map, mapping to ", $spp1_introns[$h], "\n" if $debug > 3;
            }
         }
      }
   }
}

#Now we output the orthologous intron groups:
print STDOUT "OG";
for my $spp (keys %species_found) {
   print STDOUT "\t${spp}";
}
print STDOUT "\n";
orthogroup: for my $OG (sort { substr($a, 2) <=> substr($b, 2) } keys %intron_homology) {
intron: for my $intron (sort { $a <=> $b } keys %{$intron_homology{$OG}}) {
      my $intron_line = "${OG}.intron${intron}";
species: for my $spp (keys %species_found) {
         if (exists($intron_homology{$OG}{$intron}{$spp})) {
            $intron_line .= "\t" . $txid_map{$OG}{$spp} . ".intron" . $intron_homology{$OG}{$intron}{$spp};
         } else {
            next intron;
         }
      }
      print STDOUT $intron_line, "\n";
   }
}

#Now re-iterate through the FASTA files and sort each intron into
# its appropriate file based on the homology scan above:
#We choose to open and close file handles for each FASTA record
# since this obviates any concerns about open file handle limits.
#It may, however, lead to slower performance due to fh overhead.
#for my $fasta_path (@fasta_list) {
#   my $intron_fasta;
#   unless (open($intron_fasta, "<", $fasta_path)) {
#      print STDERR "Unable to open intron FASTA ${fasta_path}.\n";
#      exit 4;
#   }
#   my $fasta_header = "";
#   my $fasta_seq = "";
#   while (my $line = <$intron_fasta>) {
#      chomp $line;
#      if ($line =~ /^>/) { #Header line
#         if ($fasta_header ne "") {
#            #Process the prior record:
#            my @header_elems = split /\s+/, substr($fasta_header, 1);
#            if ($header_elems[0] =~ /^(\S+)\.intron(\d+)$/) {
#               my $transcript_ID = $1;
#               my $intron_num = $2;
#               unless (exists($orthogroup_map{$transcript_ID})) {
#                  print STDERR "No orthogroup found for ${transcript_ID}, skipping\n" if $debug > 1;
#                  $fasta_seq = "";
#                  $fasta_header = $line;
#                  next;
#               }
#               my @txid_elems = split /_/, $header_elems[0];
#               my $short_txid = pop @txid_elems;
#               my $species = join("_", @txid_elems);
#               my $OG = $orthogroup_map{$transcript_ID};
#               if (exists($intron_homology{$OG}{$species}{$intron_num})) {
#                  my $homologous_intron = $intron_homology{$OG}{$species}{$intron_num};
#                  #Output the intron to the file:
#                  my $output_fn = "${OG}_intron${homologous_intron}_unwrapped.fasta";
#                  print STDERR $header_elems[0], " goes in ${output_fn}\n" if $debug > 2;
#                  my $output_fh;
#                  unless(open($output_fh, ">>", $output_fn)) {
#                     print STDERR "Unable to open output file ${output_fn} for appending\n";
#                  } else {
#                     print $output_fh $fasta_header, "\n";
#                     print $output_fh $fasta_seq, "\n";
#                     close($output_fh);
#                  }
#               } else {
#                  print STDERR "No homology found for ${OG} ${species} intron ${intron_num} based on header ${fasta_header}\n" if $debug > 2;
#               }
#            } else {
#               print STDERR "Unable to parse intron FASTA header for ${fasta_path} ", $header_elems[0], ", skipping\n" if $debug;
#            }
#         }
#         $fasta_seq = "";
#         $fasta_header = $line;
#      } else { #Sequence line
#         $fasta_seq .= $line;
#      }
#   }
#   #Process the last record:
#   my @header_elems = split /\s+/, substr($fasta_header, 1);
#   if ($header_elems[0] =~ /^(\S+)\.intron(\d+)$/) {
#      my $transcript_ID = $1;
#      my $intron_num = $2;
#      print STDERR "No orthogroup found for ${transcript_ID}, skipping\n" if (!exists($orthogroup_map{$transcript_ID}) and $debug > 1);
#      if (exists($orthogroup_map{$transcript_ID})) {
#         my @txid_elems = split /_/, $header_elems[0];
#         my $short_txid = pop @txid_elems;
#         my $species = join("_", @txid_elems);
#         my $OG = $orthogroup_map{$transcript_ID};
#         if (exists($intron_homology{$OG}{$species}{$intron_num})) {
#            my $homologous_intron = $intron_homology{$OG}{$species}{$intron_num};
#            #Output the intron to the file:
#            my $output_fn = "${OG}_intron${homologous_intron}_unwrapped.fasta";
#            my $output_fh;
#            unless(open($output_fh, ">>", $output_fn)) {
#               print STDERR "Unable to open output file ${output_fn} for appending\n";
#            } else {
#               print $output_fh $fasta_header, "\n";
#               print $output_fh $fasta_seq, "\n";
#               close($output_fh);
#            }
#         } else {
#            print STDERR "No homology found for ${fasta_header}\n" if $debug > 2;
#         }
#      }
#   } else {
#      print STDERR "Unable to parse intron FASTA header for ${fasta_path} ", $header_elems[0], ", skipping\n" if $debug;
#   }
#   close($intron_fasta);
#}

exit 0;
