#!/usr/bin/perl
use POSIX;
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $SCRIPTNAME = "remapOrthogroups.pl";
my $VERSION = "1.2";
#Version 1.1 fixes a bug in removing carriage returns from header line
#Version 1.2 adds some debugging output, and fixes a trailing
#            whitespace/delimiter issue

=pod

=head1 NAME

remapOrthogroups.pl - Converts OrthoFinder Orthogroups.csv to original annotation IDs

=head1 SYNOPSIS

remapOrthogroups.pl [options]

 Options:
  --help,-h,-?         Display this help documentation
  --prot_id_map,-i     Input concatenated map of OrthoFinder protein IDs to
                       original annotation/transcript IDs
                       (column 1 is OrthoFinder ID, column 2 is transcript ID)
  --orthogroups,-g     Path to Orthogroups.csv file produced by OrthoFinder
  --debug,-d           Output extra information to STDERR

=head1 DESCRIPTION

remapOrthogroups.pl converts the protein IDs in the Orthogroups.csv
output file from OrthoFinder back into transcript IDs from the
original annotations using a concatenated map produced by combining
the .map files made by prep_proteomes.mk's proteome target.

=cut


my $prot_id_map = "";
my $orthogroups = "";
my $help = 0;
my $man = 0;
my $debug = 0;
my $dispversion = 0;

GetOptions('prot_id_map|i=s' => \$prot_id_map, 'orthogroups|g=s' => \$orthogroups, 'debug|d+' => \$debug, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -output => \*STDERR, -verbose => 2) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

my $mapfh;
if (! -e ${prot_id_map}) {
   print STDERR "Unable to open protein ID map file ${prot_id_map} -- does not exist.\n";
   exit 1;
}
open($mapfh, "<", $prot_id_map) or die "Unable to open protein ID map file ${prot_id_map}";
my %prot_map = ();
while (my $line = <$mapfh>) {
   chomp $line;
   my @mapelems = split /\t/, $line;
   #Map is from OrthoFinder input ID (key) to original annotation/transcriptome ID (value):
   $prot_map{${mapelems[1]}} = ${mapelems[0]};
   print STDERR "Added ", join(" => ", ${mapelems[1]}, ${mapelems[0]}), " to OrthoFinderID=>protID map\n" if $debug > 1;
}
close($mapfh);
$/ = "\r\n";

my $orthogroupsfh;
if (! -e ${orthogroups}) {
   print STDERR "Unable to open OrthoFinder Orthogroups.csv file ${orthogroups} -- does not exist.\n";
   exit 2;
}
open($orthogroupsfh, "<", $orthogroups) or die "Unable to open OrthoFinder Orthogroups.csv file ${orthogroups}";
my $headerline = <$orthogroupsfh>;
chomp $headerline;
print $headerline, "\n";
#Keep track of the number of proteomes involved:
my @header_parts = split /\t/, $headerline;
my $num_parts = scalar(@header_parts);
print STDERR "Header indicates ${num_parts} columns\n" if $debug;
while (my $line = <$orthogroupsfh>) {
   chomp $line;
   my @speciesarr = split /\t/, $line, $num_parts; #Force $num_parts to avoid truncation of trailing null elements
   my $num_spp = scalar(@speciesarr);
   print STDERR "Record indicates ${num_spp} columns whereas header says ${num_parts} for line:\n${line}\n" if $num_spp != $num_parts and $debug;
   for (my $i = 0; $i < $num_spp; $i++) {
      next if $i==0;
      my @protarr = split /, /, $speciesarr[$i];
      my $num_prots = scalar(@protarr);
      for (my $j = 0; $j < $num_prots; $j++) {
         if (exists($prot_map{$protarr[$j]})) {
            #Replace the ID if it's found in the map (it should always be):
            print STDERR "Replaced ", $protarr[$j], " with ", $prot_map{$protarr[$j]}, "\n" if $debug > 1;
            $protarr[$j] = $prot_map{$protarr[$j]};
         } else {
            print STDERR "Unable to find element in map for ${protarr[$j]}\n";
         }
      }
      #Recompile the IDs into the comma-separated list for this species:
      $speciesarr[$i]=join(",", @protarr);
   }
   #Recompile the lists of IDs for this orthogroup and print it:
   print join("\t", @speciesarr), "\n";
}
close($orthogroupsfh);
