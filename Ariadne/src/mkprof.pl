#! /usr/bin/perl -w

use strict;
use Getopt::Std;
my %opts;
getopts('trn:o:w:c:',\%opts);

my $name = "";
$name = $opts{'n'} if defined($opts{'n'});
my $outputdir = "";
$outputdir = $opts{'o'} if defined($opts{'o'});
my $lastaidafile = "";
my $weights = "";
$weights = $opts{'w'} if defined($opts{'w'});
my $combinations = "30:100";
$combinations = $opts{'c'} if defined($opts{'c'});

while ( <> ) {
  if ( /^>> (\S+) sub-run number (\d+)/ ) {
    if ( !$outputdir ) {
      $outputdir = "$1.prof";
    }
    if ( !$name ) {
	$name = $1;
    }
    my $dir = "$outputdir/mc/$name-$2";
    my $aidafile = "$1:$2.aida";
    my $copy = "cp";
    if ( ! -f $aidafile ) {
	 $aidafile = "$1:$2.yoda";
	 $copy = "yoda2aida";
    }
    $lastaidafile = "$dir/out.aida";
    system("mkdir -p $dir");
    if ( ! -f $aidafile ) {
	print STDERR
	    "Cannot find aida/yoda file, assuming file has already been copied.";
	last;
    } else {
	system("$copy $aidafile $dir/out.aida");
    }

    if ( defined($opts{'r'}) ) {
      system("rm $aidafile");
    }
    system("rm -f $dir/used_params");
    while ( <> ) {
      if ( /^\s+set\s+(\S+)\s+(\S+)/ ) {
	system("echo $1 $2 >> $dir/used_params");
      } else {
	last;
      }
    }
  }
}



system("mkdir -p $outputdir/ref");
my $rivetdata = `rivet-config --datadir`;
chomp($rivetdata);
my %done;
open(AIDA,"<$lastaidafile") or die "No file!";
while ( <AIDA> ) {
  if ( /path="\/([^"]+)"/ ) {
      if ( !exists($done{$1}) ) {
	  my $rivetfile = "$1.aida";
	  my $copy = "cp";
	  if ( ! -f "$rivetdata/$rivetfile" ) {
	      $rivetfile = "$1.yoda";
	      $copy = "yoda2aida";
	      system("$copy $rivetdata/$rivetfile $outputdir/ref/$1.aida");
	      $done{$1} = "done";
	  }
      }
  }
}

if ( $weights ) {
  die "Cannot find weights file $weights" if system("cp $weights $outputdir/weights");
} else {
  system("cd $outputdir; prof-lsobs --mcdir mc --weight=1.0 > weights");
}
system("cd $outputdir; prof-runcombs --mcdir mc -c $combinations -o runcombs.dat; prof-interpolate --datadir . --weights weights --runs runcombs.dat");
if ( defined($opts{'t'}) ) {
  system("cd $outputdir; prof-tune --datadir . --weights weights --runs runcombs.dat");
}


