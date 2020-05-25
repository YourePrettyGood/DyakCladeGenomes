#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

#This script takes a BED6 from a reference species, where the Name (column 4)
# contains OGID (orthogroup ID), and an equivalent BED6 from a query species,
# where each line represents that species' homolog from a single-copy
# orthogroup. The BED6 Name column also should include bands, which is a
# comma-separated list of cytological bands overlapping the locus.
#From there, we iterate through the contigs of the query species, and output
# the cytological intervals spanned by subsequences of the contigs.
#These cytological intervals can then be compared to those derived from
# chromosome squashes, and help in scaffolding based on cytological evidence.

my $SCRIPTNAME = "findAlignedSCOsequences.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

findAlignedSCOsequences.pl - Find collinear sequences of SCOs between two assemblies

=head1 SYNOPSIS

findAlignedSCOsequences.pl -r [reference BED6] -q [query BED6]

 Options:
  --help,-h,-?	        Display this help documentation
  --version,-v          Output version string
  --debug,-d            Increase verbosity of output
  --query,-q            Query assembly BED6
  --reference,-r        Reference assembly BED6

=head1 DESCRIPTION

This script takes a BED6 from a reference species, where the Name (column 4)
 contains OGID (orthogroup ID), and an equivalent BED6 from a query species,
 where each line represents that species' homolog from a single-copy
 orthogroup. The BED6 Name column also should include bands, which is a
 comma-separated list of cytological bands overlapping the locus.

From there, we iterate through the contigs of the query species, and output
 the cytological intervals spanned by subsequences of the contigs.

These cytological intervals can then be compared to those derived from
 chromosome squashes, and help in scaffolding based on cytological evidence.

=cut

#Initialize the input parameters
my $help = 0;
my $man = 0;
my $debug = 0;
my $query_bed = "";
my $ref_bed = "";
my $dispversion = 0;

#Fetch the command line parameters
GetOptions('query|q=s' => \$query_bed, 'reference|r=s' => \$ref_bed, 'version|v' => \$dispversion, 'help|h|?+' => \$help, 'debug|d+' => \$debug, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

pod2usage(-exitval => 2, -output => \*STDERR) unless -e $query_bed and -e $ref_bed;

#Functions:
#Compare two cytological bands to see which comes first
sub bandDistance {
   my $a = shift @_;
   my $b = shift @_;
   my $direction = shift @_;

   #Constants used for NA values:
   my $NA_SMALL = 0;
   my $NA_LARGE = 255; #Using 255 since that's max extended ASCII value
   
   #Break the inputs apart into band ID components:
   my ($a_band, $a_subband, $a_subsubband);
   my ($b_band, $b_subband, $b_subsubband);
   if ($a =~ /^(\d+)([A-Z]?)(\d*)$/) {
      ($a_band, $a_subband, $a_subsubband) = ($1, $2, $3);
   } else {
      print STDERR "Unable to break band $a into parts, did not match regex.\n";
      return undef;
   }
   if ($b =~ /^(\d+)([A-Z]?)(\d*)$/) {
      ($b_band, $b_subband, $b_subsubband) = ($1, $2, $3);
   } else {
      print STDERR "Unable to break band $b into parts, did not match regex.\n";
      return undef;
   }

   #We address missing subbands and subsubbands by substituting with 0 or max
   # possible value depending on $direction:
   my $na_val;
   if ($direction > 0) {
      $na_val = $NA_LARGE;
   } else {
      $na_val = $NA_SMALL;
   }
   $a_subband = chr($na_val) if $a_subband eq "";
   $b_subband = chr($na_val) if $b_subband eq "";
   $a_subsubband = $na_val if $a_subsubband eq "";
   $b_subsubband = $na_val if $b_subsubband eq "";
   if ($a_band != $b_band) {
      return $a_band - $b_band;
   } else { #We need to examine further ID components
      if (ord($a_subband) != ord($b_subband)) {#ord() returns 0 for empty string, so still works in that case
         return ord($a_subband) - ord($b_subband);
      } else { #Matching up to subband (subband may be absent from both and still get here)
         return $a_subsubband - $b_subsubband;
      }
   }
}

#sub compareBands {
#   my $a = shift @_;
#   my $b = shift @_;
#
#   my $band_dist = bandDistance($a, $b);
#   if (defined($band_dist)) {
#      if ($band_dist > 0) {
#         return 1;
#      } elsif ($band_dist == 0) {
#         return 0;
#      } else {
#         return -1;
#      }
#   } else {
#      return undef;
#   }
#}

#Take an array of cytological bands, assuming no breaks in the alignment,
# and reduces redundant bands to a single instance.
#sub uniqueBands {
#   my @bands = @{shift @_};
#   my %band_hash = ();
#   my @unique_bands = ();
#
#   for my $band (@bands) {
#      push @unique_bands, $band unless exists($band_hash{$band});
#      $band_hash{$band} = undef;
#   }
#   return \@unique_bands;
#}

#Take an array of cytological bands and an orientation for the path,
# then perform a reduction on the array, eliminating redundancy and
# consolidating adjacent bands into an interval.
#sub consolidateBands {
#   my @bands = @{shift @_};
#   my $orient = shift @_;
#
#   my @unique_bands = @{uniqueBands(@bands)};
#   my $num_unique_bands = scalar(@unique_bands);
#   my @reduced_bands = ();
#   push @reduced_bands, $unique_bands[0];
#   for (my $i = 1; $i < $num_unique_bands; $i++) {
#      
#   
#   return \@reduced_bands;
#}

#Given an array of bands derived from an aligned set of SCOs, and the
# orientation of the alignment, take the most precise start and end bands
# as the consolidated interval.
sub simpleInterval {
   my @bands = @{shift @_};
   my $orient = shift @_;

   return undef unless scalar(@bands) > 0;
   my ($min_band, $max_band) = ($bands[0], $bands[0]);
   for my $band (@bands) {
#      print STDERR join("\t", $min_band, $band, bandDistance($min_band, $band, 1)), "\n";
#      print STDERR join("\t", $max_band, $band, bandDistance($max_band, $band, 0)), "\n";
      $min_band = $band if bandDistance($min_band, $band, 1) > 0;
      $max_band = $band if bandDistance($max_band, $band, 0) < 0;
   }

   my ($left_band, $right_band) = ("", "");
   if ($orient eq "+") {
      $left_band = $min_band;
      $right_band = $max_band;
   } elsif ($orient eq "-") {
      $left_band = $max_band;
      $right_band = $min_band;
   } else {
      print STDERR "Invalid alignment orientation $orient\n";
      return undef;
   }
   return join("-", $left_band, $right_band);
}

#Set up the hashes for the main algorithm:
my %ref_map = (); #Hash of hashes, first level key is OGID, second level is "left", "right", "scaf", or "bands"
my %query_map = (); #Hash of arrays of arrays, first level key is scaffold ID, array holds OGIDs in aligned segments
my @scaf_order = (); #Array of query scaffolds in the order they were read in
my %query_orients = (); #Hash of arrays, first level key is scaffold ID, array holds orientations of the aligned OG sequences
my %query_bands = (); #Hash of arrays, first level key is scaffold ID, array holds list of cytological bands in aligned segments

#Read through the reference BED6 and fill the ref_map:
print STDERR "Reading reference BED6 ${ref_bed}\n";
my $ref_fh;
unless(open($ref_fh, "<", $ref_bed)) {
   print STDERR "Error opening reference BED6 file ${ref_bed}\n";
   exit 2;
}
my ($prev_scaf, $prev_ogid) = ("", "");
while (my $ref_line = <$ref_fh>) {
   chomp $ref_line;
   my ($ogid, $bands);
   my ($scaf, $BEDstart, $BEDend, $name, $score, $strand) = split /\t/, $ref_line, 6;
   my @tags = split /;/, $name;
   for my $tag (@tags) {
      my @elems = split /=/, $tag;
      if ($elems[0] eq "OGID") {
         $ogid = $elems[1];
      } elsif ($elems[0] eq "bands") {
         $bands = $elems[1];
      }
   }
   $bands = "" unless defined($bands);
   my @bandlist = split /,/, $bands;
   if (defined($ogid)) {
      print STDERR "OGID ${ogid} found more than once in reference BED6 ${ref_bed}, overwriting in map.\n" if exists($ref_map{$ogid});
      $ref_map{$ogid} = {} unless exists($ref_map{$ogid});
      #If we're no longer contiguous with the previous (or at the very start):
      if ($scaf ne $prev_scaf or $prev_ogid eq "") {
         #left is undefined
         $ref_map{$ogid}{"left"} = undef;
      } else {
         #left is prev_ogid, prev_ogid's right is now ogid
         $ref_map{$ogid}{"left"} = $prev_ogid;
         $ref_map{$prev_ogid}{"right"} = $ogid;
      }
      $ref_map{$ogid}{"right"} = undef;
      $ref_map{$ogid}{"bands"} = defined($bands) ? \@bandlist : undef;
      $ref_map{$ogid}{"scaf"} = $scaf;
      #Move forward along the scaffold
      $prev_scaf = $scaf;
      $prev_ogid = $ogid;
   }
}
close($ref_fh);

#Read through the query BED6 and fill the query_map:
print STDERR "Reading query BED6 ${query_bed}\n";
my $query_fh;
unless(open($query_fh, "<", $query_bed)) {
   print STDERR "Error opening query BED6 file ${query_bed}\n";
   exit 3;
}

($prev_scaf, $prev_ogid) = (undef, undef);
while (my $query_line = <$query_fh>) {
   chomp $query_line;
   my ($ogid, $bands);
   my ($scaf, $BEDstart, $BEDend, $name, $score, $strand) = split /\t/, $query_line, 6;
   $prev_ogid = undef if defined($prev_scaf) and $scaf ne $prev_scaf;
   my @tags = split /;/, $name;
   for my $tag (@tags) {
      my @elems = split /=/, $tag;
      if ($elems[0] eq "OGID") {
         $ogid = $elems[1];
#      } elsif ($elems[0] eq "bands") {
#         $bands = $elems[1];
      }
   }
   push @scaf_order, $scaf unless exists($query_map{$scaf});
   $query_map{$scaf} = [] unless exists($query_map{$scaf});
   $query_orients{$scaf} = [] unless exists($query_orients{$scaf});
   $query_bands{$scaf} = [] unless exists($query_bands{$scaf});
   if (defined($prev_ogid)) { #If we're on the same scaffold and after the first OG
      #Check for adjacency of OGs in sequence, and determine orientation:
      my $orient = undef;
      if (defined($ref_map{$ogid}{"left"}) and $ref_map{$ogid}{"left"} eq $prev_ogid) {
         #Previous OG is left of current OG, so + orientation
         $orient = "+";
         my $pathnum = scalar(@{$query_map{$scaf}});
         push @{$query_map{$scaf}[$pathnum-1]}, $ogid;
         print STDERR "Overwriting orientation for path ", $pathnum, " as ", $orient, " rather than ", $query_orients{$scaf}[$pathnum-1], "\n" if defined($query_orients{$scaf}[$pathnum-1]) and $query_orients{$scaf}[$pathnum-1] ne $orient;
         $query_orients{$scaf}[$pathnum-1] = $orient;
         push @{$query_bands{$scaf}[$pathnum-1]}, @{$ref_map{$ogid}{"bands"}};
      } elsif (defined($ref_map{$ogid}{"right"}) and $ref_map{$ogid}{"right"} eq $prev_ogid) {
         #Previous OG is right of current OG, so - orientation
         $orient = "-";
         my $pathnum = scalar(@{$query_map{$scaf}});
         push @{$query_map{$scaf}[$pathnum-1]}, $ogid;
         print STDERR "Overwriting orientation for path ", $pathnum, " as ", $orient, " rather than ", $query_orients{$scaf}[$pathnum-1], "\n" if defined($query_orients{$scaf}[$pathnum-1]) and $query_orients{$scaf}[$pathnum-1] ne $orient;
         $query_orients{$scaf}[$pathnum-1] = $orient;
         push @{$query_bands{$scaf}[$pathnum-1]}, @{$ref_map{$ogid}{"bands"}};
      } else {
         #No match found, so current OG is the start of a new segment
         push @{$query_map{$scaf}}, [$ogid];
         push @{$query_bands{$scaf}}, [];
         my $pathnum = scalar(@{$query_map{$scaf}});
         push @{$query_bands{$scaf}[$pathnum-1]}, @{$ref_map{$ogid}{"bands"}};
      }
   } else {
      #First OG of the scaffold, so add it into the query map
      push @{$query_map{$scaf}}, [$ogid];
      push @{$query_bands{$scaf}}, [];
      my $pathnum = scalar(@{$query_map{$scaf}});
      push @{$query_bands{$scaf}[$pathnum-1]}, @{$ref_map{$ogid}{"bands"}};
   }
   $prev_scaf = $scaf;
   $prev_ogid = $ogid;
}
close($query_fh);

#Generate non-redundant list of cytological bands for each scaffold:
#for my $scaf (keys %query_map) {
for my $scaf (@scaf_order) {
   next unless exists($query_map{$scaf});
   my $pathnum = scalar(@{$query_map{$scaf}});
   for (my $i = 0; $i < $pathnum; $i++) {
      my $pathlength = scalar(@{$query_map{$scaf}[$i]});
      unless (defined($query_orients{$scaf}[$i])) {
         print STDERR "Single SCO interval for ${scaf} interval ${i}, unable to infer orientation.\n" if $debug;
         print STDOUT join("\t", $scaf, $i, join(",", @{$ref_map{$query_map{$scaf}[$i][0]}{"bands"}}), ".", $pathlength), "\n";
      } else {
         my $interval = simpleInterval(\@{$query_bands{$scaf}[$i]}, $query_orients{$scaf}[$i]);
         $interval = "" unless defined($interval);
         print STDOUT join("\t", $scaf, $i, $interval, $query_orients{$scaf}[$i], $pathlength), "\n";
      }
      if ($debug) {
         for my $ogid (@{$query_map{$scaf}[$i]}) {
            print STDOUT join("\t", $scaf, $ogid, @{$ref_map{$ogid}{"bands"}}), "\n"; #Debugging
         }
      }
   }
}

exit 0;
