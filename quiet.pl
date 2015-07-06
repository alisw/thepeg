#!/usr/bin/perl -pi

if ( /\@LIBTOOL\@/ ) {
  print "ifdef SHOWCOMMAND\n";
  print;
  print "else\n";
  s/\@LIBTOOL\@/\@LIBTOOL\@ --quiet/;
  print;
  print "endif\n";
  $_ = "";
}
elsif ( /^\s+\$\(CXXLINK\)/ ) {
  print "ifdef SHOWCOMMAND\n";
  print;
  print "else\n\t\@echo \"linking \$(subdir)/\$\@ ...\"\n";
  s/\$\(CXXLINK\)/\@\$\(CXXLINK\)/;
  print;
  print "endif\n";
  $_ = "";
}
elsif ( /^\s+\$\(LINK\)/ ) {
  print "ifdef SHOWCOMMAND\n";
  print;
  print "else\n\t\@echo \"linking \$(subdir)/\$\@ ...\"\n";
  s/\$\(LINK\)/\@\$\(LINK\)/;
  print;
  print "endif\n";
  $_ = "";
}
elsif ( /^\s+\$\(\w+_LINK\)/ ) {
  print "ifdef SHOWCOMMAND\n";
  print;
  print "else\n\t\@echo \"linking \$(subdir)/\$\@ ...\"\n";
  s/^\t\$/\t\@\$/;
  print;
  print "endif\n";
  $_ = "";
}
elsif ( /^\s+\$\(F77COMPILE\)/ ) {
  print "ifdef SHOWCOMMAND\n";
  print;
  print "else\n\t\@echo \"compiling \$(subdir)/\$\@ ...\"\n";
  s/\$\(F77COMPILE\)/\@\$\(F77COMPILE\)/;
  print;
  print "endif\n";
  $_ = "";
}
elsif ( /^\s+\$\(LTF77COMPILE\)/ ) {
  print "ifdef SHOWCOMMAND\n";
  print;
  print "else\n\t\@echo \"compiling \$(subdir)/\$\@ ...\"\n";
  s/\$\(LTF77COMPILE\)/\@\$\(LTF77COMPILE\)/;
  print;
  print "endif\n";
  $_ = "";
}
elsif ( /\.cc\.lo/ || /\.cc\.o/ || /\.cc\.obj/ || /^lib\S+\.lo:\s+\S+\.cc/ ||
	/^lib\S+\.lo:\s+\S+\.cpp/ ) {
  $chunk = "else\n\t\@echo \"compiling \$(subdir)/\$< ...\"\n";
  print;
  print "ifdef SHOWCOMMAND\n";
  $_ = "";
}
elsif ( $_ =~ /^\s*$/ && $chunk ) {
  $chunk .= "endif\n";
  $chunk =~ s/if \$\(LTCXXCOMPILE\)/\@if \$\(LTCXXCOMPILE\)/g;
  $chunk =~ s/if \$\(CXXCOMPILE\)/\@if \$\(CXXCOMPILE\)/g;
  $chunk =~ s/@am__fastdepCXX_TRUE@\t\$\(LTCXXCOMPILE\)/@am__fastdepCXX_TRUE@\t\@\$\(LTCXXCOMPILE\)/g;
  $chunk =~ s/@am__fastdepCXX_TRUE@\t\$\(CXXCOMPILE\)/@am__fastdepCXX_TRUE@\t\@\$\(CXXCOMPILE\)/g;
  $chunk =~ s/@am__fastdepCXX_TRUE@\t\$\(LIBTOOL\)/@am__fastdepCXX_TRUE@\t\@\$\(LIBTOOL\)/g;
  $chunk =~ s/@am__fastdepCXX_TRUE@\tif \$\(LIBTOOL\)/@am__fastdepCXX_TRUE@\t\@if \$\(LIBTOOL\)/g;
  $chunk =~ s/@am__fastdepCXX_TRUE@\tmv -f/@am__fastdepCXX_TRUE@\t\@mv -f/g;
  print $chunk;
  $chunk = "";
}
elsif ( $chunk ) {
  $chunk .= $_;
}
