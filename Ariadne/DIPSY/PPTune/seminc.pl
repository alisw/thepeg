#!/usr/bin/perl -w

use strict;
sub inputerror($);
sub writefile($$$$$);

my @energy;
my @frame;
my @xsec;
my @xsecerr;
my @elxsec;
my @elxsecerr;
my $i = 0;
my $runname = "";
my $E = 0;
my $hist = "d00";
my $x = 0;
my $dx = 0;

while ( <> ) {

    if ( />> (\S+) sub-run number\s+(\d+)/ ) {
	$i = $2;
	$runname = $1;
    }
    elsif ( /Starting DIPSY run at CoM energy (\d+) GeV/ ) {
	$E = $1;
	if ( $E == 7000 ) {
	    $hist = "d01";
	}
	elsif ( $E == 1800 ) {
	    $hist = "d02";
	}
	elsif ( $E == 540 ) {
	    $hist = "d03";
	}
	elsif ( $E == 140 ) {
	    $hist = "d04";
	}
	else {
	    $hist = "";
	}
    }
    elsif ( /SemiIncl: Total:\s*(\S+)/ ) {
	if ( $hist && $hist ne "d03" && $runname ne "" ) {
	    ($x, $dx) = inputerror($1);
	    writefile("$runname:$i.yoda", "$hist-x01-y01", $E, $x, $dx);
	}
    }
    elsif ( /SemiIncl: Elastic:\s*(\S+)/ ) {
	if ( $hist && $hist ne "d04" && $runname ne "" ) {
	    ($x, $dx) = inputerror($1);
	    writefile("$runname:$i.yoda", "$hist-x01-y02", $E, $x, $dx);
	}
    }
}

if ( $runname eq "" ) {
    exit(1);
}

sub inputerror($) {


  if ( $_[0] =~ /([-+]?\d+\.?(\d*))\((\d+)\)([eE][+-]?\d+)/ ) {
    return ( "$1$4", $3*"1$4"*10**(-length($2)) );
  }
  return (0, 0);

}


sub writefile($$$$$) {
    my $file = $_[0];
    my $tag = $_[1];
    my $E = $_[2];
    my $x = $_[3];
    my $dx = $_[4];
    my $old = 0;

    if ( -e $file ) {
	system("mv $file $file~" );
	$old = 1;
    }
    open NFL, ">$file";
    if ( $old ) {
	open OFL, "<$file~";
	my $keep = 1;
	while ( <OFL> ) {
	    if ( $keep ) {
		if ( /# BEGIN YODA_SCATTER2D \/TOTALXSEC\/$tag/ ) {
		    $keep = 0;
		} else {
		    print NFL;
		}
	    } else {
		if ( /# END YODA_SCATTER2D/ ) {
		    $keep = 1;
		    <OFL>;
		}
	    }
        }
    } else {
	open NFL, ">$file";
    }
    print NFL "# BEGIN YODA_SCATTER2D /TOTALXSEC/$tag\n";
    print NFL "Path=/TOTALXSEC/$tag\n";
    print NFL "Type=Scatter2D\n";
    print NFL "# xval   xerr-   xerr+   yval    yerr-   yerr+\n";
    printf NFL "%e    %e    %e    %e    %e    %e\n",
    $E, 5.0, 5.0, $x/1000000.0, $dx/1000000.0, $dx/1000000.0;
    print NFL "# END YODA_SCATTER2D\n\n";
    close NFL;
    if ( $old ) {
	close OFL;
	system("rm $file~" );
    }
}

