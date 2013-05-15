#!/usr/bin/perl -w

#
# Script to sync matefiles for PE alignment using Illumina data
# Agilent Technologies
#

use strict;
use warnings;
use Getopt::Long;
use FileHandle;
use IO::Uncompress::Gunzip;

my($in1,$in2,$out1,$out2);
GetOptions("i=s" => \$in1,"j=s" => \$in2,"o=s" => \$out1,"p=s" => \$out2);
usage() if(!$in1 || !$in2 || !$out1 || !$out2);


# read file
my %h;


# hash both mates
hashFile($in1,"m1");
hashFile($in2,"m2");


# define filehandles for printing mates
if ($out1 =~ /.gz$/) {
  open(MATE1, "| gzip -c > $out1") || die "Can't open file";
} else {
  open(MATE1, "> $out1") or die $!;
}
if ($out2 =~ /.gz$/) {
  open(MATE2, "| gzip -c > $out2") || die "Can't open file";
} else {
  open(MATE2, "> $out2") or die $!;
}

# loop over hash and print to each output
for my $id (keys %h) {
  if($h{$id}{'m1'}{'seq'} && $h{$id}{'m2'}{'seq'}) {
    print MATE1 printRead($h{$id}{'m1'});
    print MATE2 printRead($h{$id}{'m2'});
  }
}


close(MATE1);
close(MATE2);


# print the fastq read
sub printRead {
  my %r=%{$_[0]};
  return "$r{'id'}\n$r{'seq'}\n+\n$r{'qual'}\n";
}


# set id to hash
sub hashFile {
  my($file,$mate)=@_;
  my($id, $rid, $seq, $qual);
  if ($file =~ /.gz$/) {
    open(FH, "gunzip -c $file |") || die "Can't open file";
  } else {
    open(FH, $file) || die "Can't open file";
  }
  while(<FH>) {
    chomp;
    $id = $_;
    ($rid,undef) = $id =~ m/^(.*?)[\s\/][12].*?$/;
    $seq = <FH>; chomp($seq);
      <FH>; # Jump one row...
    $qual = <FH>; chomp($qual);
    if(length($seq) != 0) {
      $h{$rid}{$mate}{'id'}=$id;
      $h{$rid}{$mate}{'seq'}=$seq;
      $h{$rid}{$mate}{'qual'}=$qual;
    }
  }
  close(FH);
}


# usage information
sub usage {
  print "\nUsage: $0\n 
-i     Input for mate1
-j     Input for mate2
-o     Output for mate1
-p     Output for mate2\n\n";
  exit(1);
}
