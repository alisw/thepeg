use strict;

sub inputerror($) {


  if ( $_[0] =~ /([-+]?\d+\.?(\d*))\((\d+)\)([eE][+-]?\d+)/ ) {
    return ( "$1$4", $3*"1$4"*10**(-length($2)) );
  }
  return (0, 0);

}

1;

