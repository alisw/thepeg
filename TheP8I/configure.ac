dnl Process this file with autoconf to produce a configure script.
dnl

AC_PREREQ([2.59])
AC_INIT([TheP8I],[0.9.a],[http://www.thep.lu.se/TheP8I/],[TheP8I])
AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_MACRO_DIR([m4])
THEPEG_LIBTOOL_VERSION_INFO(3,0,0)

AC_CONFIG_SRCDIR([Strategy/TheP8IStrategy.h])
AC_CONFIG_HEADERS([Config/config.h])
AC_LANG([C++])
AM_INIT_AUTOMAKE([1.9 gnu dist-bzip2 -Wall])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET #???
AC_PROG_LN_S

LT_PREREQ([2.2])
LT_INIT([disable-static dlopen pic-only])

THEPEG_SEARCH_THEPEG

THEPEG_SEARCH_LHAPDF

PYTHIA8WORKING="no"

AC_ARG_WITH(pythia8,
        AC_HELP_STRING([--with-pythia8=DIR],[location of Pythia8 installation]),
        [],
	[if test -d ${prefix}/include/Pythia8; then
	    with_pythia8="${prefix}"
	 else
	    with_pythia8="${prefix}/lib/Pythia8"
	 fi])
PYTHIA8_DIR=${with_pythia8}
PYTHIA8_LDFLAGS="-Wl,-rpath -Wl,$PYTHIA8_DIR/lib -L$PYTHIA8_DIR/lib"
PYTHIA8_INCLUDE="-I$PYTHIA8_DIR/include"
if test -d $PYTHIA8_DIR/include/Pythia8; then
   PYTHIA8_INCLUDE="-I$PYTHIA8_DIR/include -I$PYTHIA8_DIR/include/Pythia8"
fi
PYTHIA8_LIBS="-lpythia8"
if test "x$HAS_LHAPDF" == "no"; then
   PYTHIA8_LIBS="$PYTHIA8_LIBS -llhapdfdummy"
else
   PYTHIA8_LIBS="$PYTHIA8_LIBS $LHAPDF_LIBS"
   LDFLAGS="$LDFLAGS $LHAPDF_LDFLAGS"
fi

AC_CHECK_FILE($PYTHIA8_DIR/share/Pythia8/xmldoc/Index.xml,
              PYTHIA8_XMLDIR="$PYTHIA8_DIR/share/Pythia8/xmldoc",
              AC_CHECK_FILE($PYTHIA8_DIR/xmldoc/Index.xml,
                             PYTHIA8_XMLDIR="$PYTHIA8_DIR/xmldoc",
                             AC_MSG_ERROR([
Cannot build TheP8I without a working Pythia8 xmldoc directory.
Please supply an argument: --with-pythia8=/path/to/Pythia8/installation.])))

LIBS="$LIBS $PYTHIA8_LIBS"
LDFLAGS="$LDFLAGS $PYTHIA8_LDFLAGS"
CPPFLAGS="$CPPFLAGS $PYTHIA8_INCLUDE -DPYTHIA8_XMLDIR=\"\\\"$PYTHIA8_XMLDIR\\\"\""

AC_MSG_CHECKING([if Pythia8 library and headers are present and working])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include "Pythia.h"]],
                                  [[Pythia8::Pythia pythia;]])],
				  AC_MSG_RESULT([$PYTHIA8_DIR]) &&
				  PYTHIA8WORKING="yes",
   				  AC_MSG_RESULT([NOT FOUND]) &&
				  AC_MSG_ERROR([
Cannot build TheP8I without a working Pythia8 installation.
Please supply an argument: --with-pythia8=/path/to/Pythia8/installation.]))


AC_PATH_PROG(PERL, perl)


AM_CXXFLAGS="-ansi  -Wall"
AC_SUBST(AM_CXXFLAGS)
AM_CONDITIONAL([PYTHIA8WORKS], [test "x$PYTHIA8WORKING" == "xyes"])
AC_SUBST(PYTHIA8_DIR)
AC_SUBST(PYTHIA8_XMLDIR)
THEPEG_EMPTY_SUBST

AC_CONFIG_FILES([Config/Makefile
include/Makefile
Strategy/Makefile
Hadronization/Makefile
BE/Makefile
lib/Makefile
src/Makefile
Doc/Makefile
Doc/refman.conf
Doc/refman.h
Makefile])

AC_CONFIG_FILES([src/P8Particles.pl],[chmod +x src/P8Particles.pl])
AC_CONFIG_FILES([Config/interfaces.pl],[chmod +x Config/interfaces.pl])
AC_CONFIG_FILES([BE/BoseEinsteinHandler-init.h],[Config/interfaces.pl --init BE/BoseEinsteinHandler BoseEinstein])
AC_CONFIG_FILES([BE/BoseEinsteinHandler-var.h],[Config/interfaces.pl --var BE/BoseEinsteinHandler BoseEinstein])
AC_CONFIG_FILES([BE/BoseEinsteinHandler-input.h],[Config/interfaces.pl --input BE/BoseEinsteinHandler BoseEinstein])
AC_CONFIG_FILES([BE/BoseEinsteinHandler-output.h],[Config/interfaces.pl --output BE/BoseEinsteinHandler BoseEinstein])
AC_CONFIG_FILES([BE/BoseEinsteinHandler-interfaces.h],[Config/interfaces.pl --interfaces BE/BoseEinsteinHandler BoseEinstein])


AC_CONFIG_FILES([Hadronization/StringFragmentation-init.h],[Config/interfaces.pl --init Hadronization/StringFragmentation StringFragmentation StringZ StringPT StringFlav FragmentationSystems MiniStringFragmentation HadronLevel])
AC_CONFIG_FILES([Hadronization/StringFragmentation-var.h],[Config/interfaces.pl --var Hadronization/StringFragmentation StringFragmentation StringZ StringPT StringFlav FragmentationSystems MiniStringFragmentation HadronLevel])
AC_CONFIG_FILES([Hadronization/StringFragmentation-input.h],[Config/interfaces.pl --input Hadronization/StringFragmentation StringFragmentation StringZ StringPT StringFlav FragmentationSystems MiniStringFragmentation HadronLevel])
AC_CONFIG_FILES([Hadronization/StringFragmentation-output.h],[Config/interfaces.pl --output Hadronization/StringFragmentation StringFragmentation StringZ StringPT StringFlav FragmentationSystems MiniStringFragmentation HadronLevel])
AC_CONFIG_FILES([Hadronization/StringFragmentation-interfaces.h],[Config/interfaces.pl --interfaces Hadronization/StringFragmentation StringFragmentation StringZ StringPT StringFlav FragmentationSystems MiniStringFragmentation HadronLevel])


AC_OUTPUT

