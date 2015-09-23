dnl Process this file with autoconf to produce a configure script.
dnl

AC_PREREQ([2.59])
AC_INIT([Pythia7],[1.7.0],[http://www.thep.lu.se/Pythia7/],[Pythia7])
AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_MACRO_DIR([m4])
THEPEG_LIBTOOL_VERSION_INFO(10,0,0)

AC_CONFIG_SRCDIR([Strategy/Pythia7Strategy.h])
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

AM_CXXFLAGS="-ansi  -Wall"
AC_SUBST(AM_CXXFLAGS)

THEPEG_EMPTY_SUBST
THEPEG_CHECK_RIVET

AC_CONFIG_FILES([Config/Makefile
include/Makefile
Strategy/Makefile
StringFrag/Makefile
Shower/Makefile
lib/Makefile
src/Makefile
Doc/Makefile
Doc/refman.conf
Doc/refman.h
Makefile])
## don't know what to do here  'Doc/Makefile'

AC_OUTPUT

