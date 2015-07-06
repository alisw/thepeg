#!/usr/bin/perl -w

use strict;

while ( @ARGV ) {

  print "** $ARGV[0] **\n";

  open(PIPE, "prof-showminresults $ARGV[0]|")
    or die "prof-showminresults failed";

  my $file;
  if ( $ARGV[0] =~ /(\w+).prof\// ) {
      $file = "$1.prin";
      open(FILE, ">>$file") or die "Could not open $file!";
      print FILE "# Tuning from $ARGV[0]\n"
  }


  my $sum = 0.0;
  my $sum2 = 0.0;
  my $nrun = 0;
  my $finish = 0;

  while ( <PIPE> ) {

    if ( /Goodness of fit\/Ndf:\s+(\S+)/ ) {
      $sum += $1;
      $sum2 += $1*$1;
      ++$nrun;
    }
    elsif ( /Summary of \d+ minimization results:/ ) {
      $finish = 1;
    }
    elsif ( $finish && /\S+/ ) {
      print;
      my ( $par, $val ) = split;
      print FILE "set $par $val\n" if $file;
    }
  }

  if ( $nrun > 0 ) {
    $sum /= $nrun;
    $sum2 /= $nrun;
    print "Chi2/Ndf: ", $sum, " (", sqrt($sum2 - $sum*$sum), ")\n";
  }
  print "-- $ARGV[0] --\n\n";
  shift @ARGV;

}
