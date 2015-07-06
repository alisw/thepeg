#!/usr/bin/perl

sub checkrun($$);

use strict;
use Getopt::Std;
my %opts;
getopts('p:',\%opts);

my $pattern = $opts{'p'} if defined($opts{'p'}); 

while (<>) {
  if ( /^saverun\s+(\S+)\s+(\S+)/ ) {
    my $saverun = checkrun($1, $pattern);
    print "$saverun ${1} $2\n";
  }
  elsif ( /^SAVERUN\s+(\S+)\s+(\S+)/ ) {
    my $saverun = checkrun($1, $pattern);
    print "set $2:EventHandler:YFrametest 0.5\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 140\n";
    print "$saverun ${1}01 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 540\n";
    print "$saverun ${1}05 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 1800\n";
    print "$saverun ${1}18 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
    print "$saverun ${1}70 $2\n";
    print "set $2:EventHandler:YFrametest 0.8\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 140\n";
    print "$saverun ${1}01Y8 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 540\n";
    print "$saverun ${1}05Y8 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 1800\n";
    print "$saverun ${1}18Y8 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
    print "$saverun ${1}70Y8 $2\n";
  }
  elsif ( /^SAVERUNYR\s+(\S+)\s+(\S+)/ ) {
    my $saverun = checkrun($1, $pattern);
    print "set $2:EventHandler:YFrametest -1.2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 140\n";
    print "$saverun ${1}01 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 540\n";
    print "$saverun ${1}05 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 1800\n";
    print "$saverun ${1}18 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
    print "$saverun ${1}70 $2\n";
  }
  elsif ( /^SAVERUNYM\s+(\S+)\s+(\S+)/ ) {
    my $saverun = checkrun($1, $pattern);
    print "set $2:EventHandler:YFrametest 0.5\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 140\n";
    print "$saverun ${1}01 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 540\n";
    print "$saverun ${1}05 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 1800\n";
    print "$saverun ${1}18 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
    print "$saverun ${1}70 $2\n";
  }
  elsif ( /^SAVERUNFS\s+(\S+)\s+(\S+)/ ) {
    my $saverun = checkrun($1, $pattern);
    print "set $2:EventHandler:LuminosityFunction:Energy 200\n";
    print "$saverun ${1}02 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 900\n";
    print "$saverun ${1}09 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 2760\n";
    print "$saverun ${1}28 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
    print "$saverun ${1}71 $2\n";
  }
  elsif ( /^SAVERUNTUNE\s+(\S+)\s+(\S+)\s+(\S+)/ ) {
      print "# expanded from SAVERUNTUNE $1 $2 $3\n";
    my $saverun = checkrun($1, $pattern);
    if ( -e "${3}00.prin" ) {
      print "read ${3}00.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 200\n";
      print "$saverun ${1}02 $2\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 900\n";
      print "$saverun ${1}09 $2\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 2760\n";
      print "$saverun ${1}28 $2\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
      print "$saverun ${1}71 $2\n";
    }
    if ( -e "${3}02.prin" ) {
      print "read ${3}02.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 200\n";
      print "$saverun ${1}0202 $2\n";
    }
    if ( -e "${3}09.prin" ) {
      print "read ${3}09.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 900\n";
      print "$saverun ${1}0909 $2\n";
    }
    if ( -e "${3}28.prin" ) {
      print "read ${3}28.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 2760\n";
      print "$saverun ${1}2828 $2\n";
    }
    if ( -e "${3}70.prin" ) {
      print "read ${3}70.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
      print "$saverun ${1}7170 $2\n";
    }
  }
  else {
      print;
  }

}

sub checkrun($$) {
    my $runname = $_[0];
    my $pattern = $_[1];
    if ( !$pattern ) {
	return "saverun";
    }
    if ( $pattern eq "all" ) {
	return "saverun";
    }
    if ( $runname =~ /$pattern/ ) {
	return "saverun";
    }
    return "# saverun";
}


