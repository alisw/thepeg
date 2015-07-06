#! /usr/bin/perl

sub inputerror($);

$analyses ="(Glauber|ElXSec|TotXSec|SemiIncl)";

@tags;
%n;
%space;
%sum;
%sum2;
%sumerr2;
%junk;

while ( <> ) {
    if ( /^($analyses:.*:)(\s*)(\S+)(.*)/ ) {
	my $tag = $1;
	my $rest = $5;
	if ( !exists($n{$tag}) ) {
	    push @tags, $tag;
	    $n{$tag} = 0;
	    $sum{$tag} = 0;
	    $sum2{$tag} = 0;
	    $sumerr2{$tag} = 0;
	    $space{$tag} = $3;
	    $junk{$tag} = $rest;
	}
	++$n{$tag};
	($x, $dx) = inputerror($4);
	$sum{$tag} += $x;
	$sum2{$tag} += $x*$x;
	$sumerr2{$tag} += $dx*$dx;
	if ( $tag =~ s/Total:/Inelastic(tot):/ ) {
	    if ( !exists($n{$tag}) ) {
		push @tags, $tag;
		$n{$tag} = 0;
		$sum{$tag} = 0;
		$sum2{$tag} = 0;
		$sumerr2{$tag} = 0;
		$space{$tag} = "      ";
	    }
	    ++$n{$tag};
	    $sum{$tag} += $x;
	    $sum2{$tag} += $x*$x;
	    $sumerr2{$tag} += $dx*$dx;
	    $junk{$tag} = $rest;
	}
	if ( $tag =~ s/Elastic:/Inelastic(tot):/ ) {
	    if ( !exists($n{$tag}) ) {
		push @tags, $tag;
		$n{$tag} = 0;
		$sum{$tag} = 0;
		$sum2{$tag} = 0;
		$sumerr2{$tag} = 0;
	    }
	    ++$n{$tag};
	    $sum{$tag} -= $x;
	    $sum2{$tag} += $x*$x;
	    $sumerr2{$tag} += $dx*$dx;
	    $junk{$tag} = $rest;
	}
    }
}

foreach $tag (@tags) {
    next if $tag =~ /Inelastic\(tot\)/;
  $x = $sum{$tag}/$n{$tag};
  $dx = sqrt($sumerr2{$tag})/$n{$tag};
  $oerr = `/home/beckett/leif/bin/outputerr $x $dx`;
  printf "%-50s%s%s\n", $tag, $oerr, $junk{$tag};
}

foreach $tag (@tags) {
    next if not $tag =~ /Inelastic\(tot\)/;
  $x = $sum{$tag}/($n{$tag}/2.0);
  $dx = sqrt($sumerr2{$tag})/($n{$tag}/2.0);
  $oerr = `/home/beckett/leif/bin/outputerr $x $dx`;
  printf "%-50s%s%s\n", $tag, $oerr, $junk{$tag};
}








sub inputerror($) {


  if ( $_[0] =~ /([-+]?\d+\.?(\d*))\((\d+)\)([eE][+-]?\d+)/ ) {
    return ( "$1$4", $3*"1$4"*10**(-length($2)) );
  }
  return (0, 0);

}
