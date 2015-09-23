#! /usr/bin/perl -w

use strict;
use Getopt::Std;

print "@ARGV\n";
my %opts;
getopts('trw:c:o:s',\%opts);

my $weights = "";
$weights = $opts{'w'} if defined($opts{'w'});
my $combinations = "30:100";
$combinations = $opts{'c'} if defined($opts{'c'});

die "No filetag specified" if !@ARGV;

my $tag = $ARGV[0];
my $profdir = "$tag.prof";
$profdir = $opts{'o'} if defined($opts{'o'});
mkdir $profdir;
mkdir "$profdir/mc";


my @files = glob("$tag\[0-9\]*.log");

foreach my $file ( @files ) {

    if ( system("../seminc.pl $file") ) {
	next;
    }
    my $YF = 0;
    if ( $file =~ /Y(\d)\.log/ ) {
	$YF = $1;
	print "Processing  $file (shifted frame $YF)... ";
    } else {
	print "Processing  $file... ";
    }


    open(LOGFILE, "<$file");

    my $nruns = 0;

    while ( <LOGFILE> ) {
	if ( /^>> (\S+) sub-run number (\d+)/ ) {
	    ++$nruns;
	    my $dir = "$profdir/mc/$tag-$2";
	    my $runno = $2;
	    my $runname = $1;
	    my $fullrunname = $1;
	    if ( $runname =~ /^(.*)#/ ) {
		$runname = $1;
	    }
	    mkdir $dir;
	    my $yodafile = "$fullrunname:$runno.yoda";
	    my $destyodafile = "$runname:$runno.yoda";
	    my $copy = "cp $yodafile $dir/$yodafile";
	    if ( $YF ) {
		system("sed -e 's/\-x0/-x8/' $yodafile > $dir/$destyodafile");
	    } else {
		system("cp $yodafile $dir/$destyodafile");
	    }
	    if ( defined($opts{'r'}) ) {
		system("rm $yodafile");
	    }
	    my $params = "$dir/used_params";
	    $params = "$dir/new_params" if  -e $params;
	    while ( <LOGFILE> ) {
		if ( /^\s+set\s+(\S+)\s+(\S+)/ ) {
		    system("echo $1 $2 >> $params");
		} else {
		    last;
		}
	    }
	    if ( -e "$dir/new_params" ) {
		if ( system("diff -q $dir/new_params $dir/used_params") ) {
		    die "Warning: inconsistent parameters in $dir/used_params";
		}
		system("rm $dir/new_params");
	    }
	}
    }
    print "$nruns runs\n";
}

foreach my $dir ( glob("$profdir/mc/$tag-*") ) {
    system("cd $dir; rm -f out.yoda; cat *.yoda > out.yoda; yoda2aida out.yoda out.aida");
}

system ("cd $profdir; rm -f ref; ln -s ../../refs ref");


if ( $weights ) {
    die "Cannot find weights file $weights"
	if system("cp $weights $profdir/weights");
}
elsif ( -e "$tag.weights" ) {
    system("cp $tag.weights $profdir/weights");
}
else {
    system("cd $profdir; prof-lsobs --mcdir mc --weight=1.0 > weights");
}

system("cd $profdir; prof-runcombs --mcdir mc -c $combinations -o runcombs.dat; prof-interpolate --datadir . --weights weights --runs runcombs.dat");
if ( defined($opts{'t'}) ) {
  system("cd $profdir; prof-tune --datadir . --weights weights --runs runcombs.dat");
}


if ( defined($opts{'s'}) ) {
    system("../profsum.pl $profdir/tunes/results.pkl");
}
