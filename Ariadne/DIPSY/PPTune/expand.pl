#!/usr/bin/perl

sub checkrun($$);
sub skiprun($$);

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
    next if skipkrun($1, $pattern);
    print "set $2:EventHandler:YFrametest 0.5\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 140\n";
    print "saverun ${1}01 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 540\n";
    print "saverun ${1}05 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 1800\n";
    print "saverun ${1}18 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
    print "saverun ${1}7O $2\n";
    print "set $2:EventHandler:YFrametest 0.8\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 140\n";
    print "saverun ${1}01Y8 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 540\n";
    print "saverun ${1}05Y8 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 1800\n";
    print "saverun ${1}18Y8 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
    print "saverun ${1}7OY8 $2\n";
  }
  elsif ( /^SAVERUNYR\s+(\S+)\s+(\S+)/ ) {
    next if skiprun($1, $pattern);
    print "set $2:EventHandler:YFrametest -1.2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 140\n";
    print "saverun ${1}01 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 540\n";
    print "saverun ${1}05 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 1800\n";
    print "saverun ${1}18 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
    print "saverun ${1}7O $2\n";
  }
  elsif ( /^SAVERUNYM\s+(\S+)\s+(\S+)/ ) {
    next if skiprun($1, $pattern);
    print "set $2:EventHandler:YFrametest 0.5\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 140\n";
    print "saverun ${1}01 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 540\n";
    print "saverun ${1}05 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 1800\n";
    print "saverun ${1}18 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
    print "saverun ${1}7O $2\n";
  }
  elsif ( /^SAVERUNY8\s+(\S+)\s+(\S+)/ ) {
    next if skiprun($1, $pattern);
    print "set $2:EventHandler:YFrametest 0.\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 140\n";
    print "saverun ${1}01 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 540\n";
    print "saverun ${1}05 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 1800\n";
    print "saverun ${1}18 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
    print "saverun ${1}7O $2\n";
    print "set $2:EventHandler:YFrametest 0.5\n";
  }
  elsif ( /^SAVERUNFS\s+(\S+)\s+(\S+)/ ) {
    next if skiprun($1, $pattern);
    print "set $2:EventHandler:LuminosityFunction:Energy 200\n";
    print "saverun ${1}02 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 900\n";
    print "saverun ${1}09 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 2760\n";
    print "saverun ${1}28 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
    print "saverun ${1}70 $2\n";
  }
  elsif ( /^SAVERUNFSX\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ ) {
    next if skiprun($1, $pattern);
    print "set $2:DefaultObjects[0] $5\n";
    print "set $4:DefaultObjects[0] $5\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 200\n";
    print "saverun ${1}02 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 900\n";
    print "saverun ${1}09 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 2760\n";
    print "saverun ${1}28 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
    print "saverun ${1}70 $2\n";
    if ( $1 eq $3 ) {
	print "set $4:EventHandler:LuminosityFunction:Energy 140\n";
	print "saverun ${1}01 $4\n";
	print "set $4:EventHandler:LuminosityFunction:Energy 540\n";
	print "saverun ${1}05 $4\n";
	print "set $4:EventHandler:LuminosityFunction:Energy 1800\n";
	print "saverun ${1}18 $4\n";
	print "set $4:EventHandler:LuminosityFunction:Energy 7000\n";
	print "saverun ${1}7O $4\n";
    } else {
	if ( -e "${3}01.log" ) {
	    system("ln -sf ${3}01.log ${1}01.log") if  ! -e "${1}01.log";
	} else {
	    print "set $4:EventHandler:LuminosityFunction:Energy 140\n";
	    print "saverun ${1}01 $4\n";
	    print "set $4:EventHandler:LuminosityFunction:Energy 7000\n";
	}
	if ( -e "${3}05.log" ) {
	    system("ln -sf ${3}05.log ${1}05.log") if  ! -e "${1}05.log";
	} else {
	    print "set $4:EventHandler:LuminosityFunction:Energy 540\n";
	    print "saverun ${1}05 $4\n";
	    print "set $4:EventHandler:LuminosityFunction:Energy 7000\n";
	}
	if ( -e "${3}18.log" ) {
	    system("ln -sf ${3}18.log ${1}18.log") if  ! -e "${1}18.log";
	} else {
	    print "set $4:EventHandler:LuminosityFunction:Energy 1800\n";
	    print "saverun ${1}18 $4\n";
	    print "set $4:EventHandler:LuminosityFunction:Energy 7000\n";
	}
	if ( -e "${3}7O.log" ) {
	    system("ln -sf ${3}7O.log ${1}7O.log") if  ! -e "${1}7O.log";
	} else {
	    print "set $4:EventHandler:LuminosityFunction:Energy 7000\n";
	    print "saverun ${1}7O $4\n";
	    print "set $4:EventHandler:LuminosityFunction:Energy 7000\n";
	}
    }
  }
  elsif ( /^SAVERUNFSXD\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ ) {
    next if skiprun("${1}XF", $pattern);
    system("mkdir -p $1");
    print "# Expanded from SAVERUNFSXD ${1}\n";
    print "set $2:DefaultObjects[0] $1\n";
    print "set $4:DefaultObjects[0] $1\n";
    print "set $2:Path $1\n";
    print "set $4:Path $1\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 200\n";
    print "saverun PPXF02 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 900\n";
    print "saverun PPXF09 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 2760\n";
    print "saverun PPXF28 $2\n";
    print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
    print "saverun PPXF70 $2\n";
    if ( $1 eq $3 ) {
	print "set $4:EventHandler:LuminosityFunction:Energy 140\n";
	print "saverun PPXF01 $4\n";
	print "set $4:EventHandler:LuminosityFunction:Energy 540\n";
	print "saverun PPXF05 $4\n";
	print "set $4:EventHandler:LuminosityFunction:Energy 1800\n";
	print "saverun PPXF18 $4\n";
	print "set $4:EventHandler:LuminosityFunction:Energy 7000\n";
	print "saverun PPXF7O $4\n";
    } else {
	if ( -e "${3}/PPXF01.log" ) {
	    system("cd $1; ln -sf ../${3}/PPXF01.log PPXF01.log")
		if  ! -e "${1}/PPXF01.log";
	} else {
	    print "set $4:EventHandler:LuminosityFunction:Energy 140\n";
	    print "saverun PPXF01 $4\n";
	    print "set $4:EventHandler:LuminosityFunction:Energy 7000\n";
	}
	if ( -e "${3}/PPXF05.log" ) {
	    system("cd $1; ln -sf ../${3}/PPXF05.log PPXF05.log")
		if  ! -e "${1}/PPXF05.log";
	} else {
	    print "set $4:EventHandler:LuminosityFunction:Energy 540\n";
	    print "saverun PPXF05 $4\n";
	    print "set $4:EventHandler:LuminosityFunction:Energy 7000\n";
	}
	if ( -e "${3}/PPXF18.log" ) {
	    system("cd $1; ln -sf ../${3}/PPXF18.log PPXF18.log")
		if  ! -e "${1}/PPXF18.log";
	} else {
	    print "set $4:EventHandler:LuminosityFunction:Energy 1800\n";
	    print "saverun PPXF18 $4\n";
	    print "set $4:EventHandler:LuminosityFunction:Energy 7000\n";
	}
	if ( -e "${3}/PPXF7O.log" ) {
	    system("cd $1; ln -sf ../${3}/PPXF7O.log PPXF7O.log")
		if  ! -e "${1}/PPXF7O.log";
	} else {
	    print "set $4:EventHandler:LuminosityFunction:Energy 7000\n";
	    print "saverun PPXF7O $4\n";
	    print "set $4:EventHandler:LuminosityFunction:Energy 7000\n";
	}
    }
  }
  elsif ( /^SAVERUNTUNE\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ ) {
    next if skiprun($1, $pattern);
    print "# expanded from SAVERUNTUNE $1 $2 $3\n";
    print "set $2:DefaultObjects[0] $4\n";
    if ( -e "${3}00.prin" ) {
      print "read ${3}00.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 200\n";
      print "saverun ${1}02 $2\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 900\n";
      print "saverun ${1}09 $2\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 2760\n";
      print "saverun ${1}28 $2\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
      print "saverun ${1}70 $2\n";
    }
    if ( -e "${3}02.prin" ) {
      print "read ${3}02.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 200\n";
      print "saverun ${1}0202 $2\n";
    }
    if ( -e "${3}09.prin" ) {
      print "read ${3}09.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 900\n";
      print "saverun ${1}0909 $2\n";
    }
    if ( -e "${3}28.prin" ) {
      print "read ${3}28.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 2760\n";
      print "saverun ${1}2828 $2\n";
    }
    if ( -e "${3}70.prin" ) {
      print "read ${3}70.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
      print "saverun ${1}7070 $2\n";
    }
  }
  elsif ( /^SAVERUNTUNED\s+(\S+)\s+(\S+)/ ) {
    next if skiprun("${1}T", $pattern);
    print "# expanded from SAVERUNTUNED $1 $2\n";
    print "set $2:DefaultObjects[0] $1\n";
    print "set $2:Path $1\n";
    if ( -e "${1}/Tune00.prin" ) {
      print "read ${1}/Tune00.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 200\n";
      print "saverun PPT02 $2\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 900\n";
      print "saverun PPT09 $2\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 2760\n";
      print "saverun PPT28 $2\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
      print "saverun PPT70 $2\n";
    }
    if ( -e "${1}/Tune02.prin" ) {
      print "read ${1}/Tune02.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 200\n";
      print "saverun PPT0202 $2\n";
    }
    if ( -e "${1}/Tune09.prin" ) {
      print "read ${1}/Tune09.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 900\n";
      print "saverun PPT0909 $2\n";
    }
    if ( -e "${1}/Tune28.prin" ) {
      print "read ${1}/Tune28.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 2760\n";
      print "saverun PPT2828 $2\n";
    }
    if ( -e "${1}/Tune70.prin" ) {
      print "read ${1}/Tune70.prin\n";
      print "set $2:EventHandler:LuminosityFunction:Energy 7000\n";
      print "saverun PPT7070 $2\n";
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

sub skiprun($$) {
    my $runname = $_[0];
    my $pattern = $_[1];
    if ( !$pattern ) {
	return 0;
    }
    if ( $pattern eq "all" ) {
	return 0;
    }
    if ( $runname =~ /$pattern/ ) {
	return 0;
    }
    return 1;
}


