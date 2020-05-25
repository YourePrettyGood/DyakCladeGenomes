#!/usr/bin/env perl
use strict;
use warnings;

my $bed_path = shift @ARGV;
my $fasta_path = shift @ARGV;

my $bedfh;
my $fastafh;
unless (open($bedfh, "<", $bed_path)) {
   print STDERR "Cannot open BED ${bed_path}\n";
   exit 1;
}
unless (open($fastafh, "<", $fasta_path)) {
   print STDERR "Cannot open FASTA ${fasta_path}\n";
   exit 2;
}

my %unmask_regions = ();
while (my $bedline = <$bedfh>) {
   my @bed_parts = split /\t/, $bedline;
   my ($scaf, $start, $end) = @bed_parts[0..2];
   $unmask_regions{$scaf} = [] unless exists($unmask_regions{$scaf});
   push @{$unmask_regions{$scaf}}, {"start"=>$start, "length"=>$end-$start};
}
close($bedfh);

my $scaffold = "";
while (my $fastaline = <$fastafh>) {
   chomp $fastaline;
   if ($fastaline =~ /^>/) {
      $scaffold = substr($fastaline, 1);
      print STDOUT $fastaline, "\n";
   } else {
      for my $regionref (@{$unmask_regions{$scaffold}}) {
         substr($fastaline, $regionref->{"start"}, $regionref->{"length"}) = uc(substr($fastaline, $regionref->{"start"}, $regionref->{"length"}));
      }
      print STDOUT $fastaline, "\n";
   }
}
close($fastafh);
