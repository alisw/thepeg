#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my $path = $0;
$path = "./$path" unless $path =~ /\//;
$path =~ s/\/[^\/]+$//;
my $doquiet = "$path/quiet.pl";

my %opts;

getopts('vcCDt:', \%opts);

my $quiet = !defined $opts{'v'};
my $doconf = defined $opts{'c'};
my $docheck = defined $opts{'C'};
my $dodistcheck = defined $opts{'D'};
if ( $docheck || $dodistcheck ) {
  $doconf = 1;
}

my $prefix = "";
$prefix = "--prefix=$ENV{'HOME'}/MCEG/Projects/Work/local "
    if -d "$ENV{'HOME'}/MCEG/Projects/Work/local";

my $tee = "";
if ( defined $opts{'t'} ) {
  system("echo \"make check\" > $opts{'t'}") == 0 or die "aborted"
    if $docheck;
  system("echo \"make distcheck\" > $opts{'t'}") == 0 or die "aborted"
    if $dodistcheck;
  $tee = "2>&1 | tee -a " . $opts{'t'} if defined $opts{'t'};
}

my @dirs;
if ( !@ARGV ) {
  my $subdirs = `grep SUBDIRS ./Makefile.am`;
  $subdirs =~ s/SUBDIRS\s*=(.*?)$/$1/;
  @dirs = split ' ', $subdirs;
} else {
  @dirs = @ARGV;
}

foreach my $dir ( @dirs ) {
  if ( "$dir" ne "ThePEG" ) {
    system("test $dir/../ThePEG/acinclude.m4 -ef $dir/acinclude.m4 || cp -u $dir/../ThePEG/acinclude.m4 $dir/acinclude.m4");
  }
  print "Configuring $dir (aclocal, ";
  system("cd $dir; aclocal --force") == 0 or die "aborted";
  print "libtoolize, ";
  system("cd $dir; libtoolize --automake") == 0 or die "aborted";
  if ( !system("grep -q AC_CONFIG_HEADERS $dir/configure.ac") ) {
    print "autoheader, ";
    system("cd $dir; autoheader") == 0 or die "aborted";
  }
  print "automake, ";
  system("cd $dir; automake  --foreign --add-missing") == 0 or die "aborted";
  print "autoconf";
  system("cd $dir; autoconf") == 0 or die "aborted";
  print ")\n";
  if ( $quiet ) {
    my $files = `find $dir -name Makefile.in`;
    $files =~ s/\n/ /gs;
    system("$doquiet $files") == 0 or die "aborted";
  }
  if ( @ARGV && $doconf ) {
    system("cd $dir; ./configure --enable-unitchecks $prefix $tee") == 0 or
	die "aborted";
    if ( $docheck ) {
      system("cd $dir; make check $tee") == 0 or die "aborted";
    }
    if ( $dodistcheck ) {
      system("cd $dir; make distcheck $tee") == 0 or die "aborted";
    }
  }
}

if ( !@ARGV ) {
  system("aclocal --force; autoconf; libtoolize --automake; automake  --foreign --add-missing") == 0 or die "aborted";
  if ( $quiet && -f "Makefile.in" ) {
    system("$doquiet Makefile.in") == 0 or die "aborted";
  }
  if ( $doconf ) {
    system("./configure --enable-unitchecks $prefix $tee") == 0 or die "aborted";
  }
  if ( $docheck ) {
    system("(time make check) $tee") == 0 or die "aborted";
  }
  if ( $dodistcheck ) {
    system("(time make distcheck) $tee") == 0 or die "aborted";
  }
}

