#!/usr/bin/perl -w

use strict;
sub inputerror($);

my @energy;
my @frame;
my @xsec;
my @xsecerr;
my @elxsec;
my @elxsecerr;
my $i;

while ( <> ) {

  if ( /sub-run number\s+(\d+)/ ) {
    $i = $1 - 1;
    push @energy, "";
    push @frame, "";
    push @xsec, "";
    push @xsecerr, "";
  }
  elsif ( /set LHCEventHandler:YFrametest (.*)/ ) {
    $frame[$i] = $1;
  }
  elsif ( /set LHCLumi:Energy (.*)/ ) {
    $energy[$i] = $1;
  }
  elsif ( /^ElXSec: Elastic/ ) {
    my @line = split;
    ($elxsec[$i], $elxsecerr[$i]) = inputerror($line[4]);
  }
  elsif ( /^TotXSec/ ) {
    my @line = split;
    ($xsec[$i], $xsecerr[$i]) = inputerror($line[4]);
  }
}

for ( $i = 0; $i < @xsec; $i += 2 ) {
  my $R = $xsec[$i]/$xsec[$i+1];
  my $dx1 = $xsecerr[$i]/$xsec[$i];
  my $dx2 = $xsecerr[$i + 1]/$xsec[$i + 1];
  my $dR = $R*sqrt($dx1*$dx1 + $dx2*$dx2);

  my $elR = $elxsec[$i]/$elxsec[$i+1];
  my $eldx1 = $elxsecerr[$i]/$elxsec[$i];
  my $eldx2 = $elxsecerr[$i + 1]/$elxsec[$i + 1];
  my $delR = $elR*sqrt($eldx1*$eldx1 + $eldx2*$eldx2);
  printf "%6s%7.1f%7.1f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f\n",
    $energy[$i], $xsec[$i]/1000000.0, $elxsec[$i]/1000000.0,
      $R, $dR, $elR, $delR, $elxsec[$i]/$xsec[$i], $xsec[$i]/$xsec[0], $elxsec[$i]/$elxsec[0] ;
}


sub inputerror($) {


  if ( $_[0] =~ /([-+]?\d+\.?(\d*))\((\d+)\)([eE][+-]?\d+)/ ) {
    return ( "$1$4", $3*"1$4"*10**(-length($2)) );
  }
  return (0, 0);

}



