#!/usr/bin/perl -w

sub inputerror($);

use strict;
use Getopt::Std;

my %opts;

getopts('rg:x:R:', \%opts);

my $ratio = defined($opts{'r'});

my $gtype = 0;
$gtype = $opts{'g'} if defined($opts{'g'});

my $xtype = "t";
$xtype = $opts{'x'} if defined($opts{'x'});

my $RATIO = $opts{'R'} if defined($opts{'R'});

my $strategy = "DIPSY";
if ( $gtype eq "b" ) {
  $strategy = "black disc";
}
elsif ( $gtype eq "g" ) {
  $strategy = "grey disc";
}
elsif ( $gtype eq "G" ) {
  $strategy = "grey3 disc";
}

my $type = "Total";
if ( $xtype eq "i" ) {
  $type = "Inelastic\\(ND\\)";
}
elsif ( $xtype eq "e" ) {
  $type = "Elastic";
}
elsif ( $xtype eq "d" ) {
  $type = "Diff. total";
}

die "No filename supplied" if !@ARGV;

foreach  my $energy ( "01", "02", "05", "10", "25" ) {

  my $filename = $ARGV[0];
  $filename =~ s/##/$energy/;

#  print "$filename\n";

#  print STDERR "Searching for Type ($strategy) in $filename\n";
  open(FILE, "<$filename") or next;


  my $x0;
  my $dx0;
  my $xi;
  my $dxi;

  while ( <FILE> ) {
    if ( /Glauber: $type: \(DIPSY\):\s*(\S+)/ ) {
      ($x0, $dx0) = inputerror($1);
    }
    if ( /Glauber: $type: \($strategy\):\s*(\S+)/ ) {
      ($xi, $dxi) = inputerror($1);
    }
  }

  close(FILE);

  if ( $RATIO ) {
    $filename = $RATIO;
    $filename =~ s/##/$energy/;
#    print STDERR "Searching for Type ($strategy) in $filename\n";
    open(FILE, "<$filename") or last;
    while ( <FILE> ) {
      if ( /Glauber: $type: \($strategy\):\s*(\S+)/ ) {
	($x0, $dx0) = inputerror($1);
      }
    }
    close(FILE);
  }

  last if !defined($xi);
  if ( $ratio or $RATIO ) {
    last if !defined($x0);
    $dxi = sqrt(($dxi/$xi)*($dxi/$xi)+($dx0/$x0)*($dx0/$x0))*$xi/$x0;
    $xi /= $x0;
  }

  print $energy*100, " $xi $dxi\n";
}

sub inputerror($) {


  if ( $_[0] =~ /([-+]?\d+\.?(\d*))\((\d+)\)([eE][+-]?\d+)/ ) {
    return ( "$1$4", $3*"1$4"*10**(-length($2)) );
  }
  return (0, 0);

}
