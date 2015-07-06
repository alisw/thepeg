dnl Process this file with autoconf to produce a configure script.
dnl

AC_PREREQ([2.59])
AC_INIT(Dummy,0.0)
AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_MACRO_DIR([m4])
THEPEG_LIBTOOL_VERSION_INFO(3,0,0)
AM_INIT_AUTOMAKE([1.9 gnu dist-bzip2 -Wall])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_LANG(C++)

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET #???
AC_PROG_LN_S

LT_PREREQ([2.2])
LT_INIT([disable-static dlopen pic-only])

THEPEG_SEARCH_LHAPDF

THEPEG_SEARCH_HEPMC

THEPEG_SEARCH_THEPEG

THEPEG_ADD_THEPEG_PATH(PYTHIA7LIBPATH,Pythia7)

AC_PROG_FC
AC_FC_LIBRARY_LDFLAGS
AC_MSG_CHECKING([if HZTool is present and works])
if test "x$HZTOOLDIR" == "x"; then
  HZTOOLDIR="/usr/local/lib"
fi
if test "x$CERNLIBDIR" == "x"; then
  if test -d /usr/local/cern/lib; then
    CERNLIBDIR="/usr/local/cern/lib"
  else
    CERNLIBDIR="/usr/local/lib"
  fi
fi
oldLIB="$LIBS"
oldLDFLAGS="$LDFLAGS"
LIBS="$LIBS -lhztool -lkernlib -lmathlib -lpacklib $FCLIBS"		
LDFLAGS="$LDFLAGS -L$HZTOOLDIR -L$CERNLIBDIR"
AC_LINK_IFELSE([AC_LANG_PROGRAM([[extern "C" { void hzh9807018_(const int&); void structm_() {} void pdfset_() {} void lulist_() {} void pylist_() {} void luexec_() {} void kpisflag_() {} }]],
                                  [[hzh9807018_(1);]])], HAS_HZTOOL="yes",
                                  HAS_HZTOOL="no")
if test "$HAS_HZTOOL" == "yes"; then
  HZTOOL_LIBS="-lhztool -lkernlib -lmathlib -lpacklib $FCLIBS"
  HZTOOL_LDFLAGS=" -L$HZTOOLDIR -L$CERNLIBDIR"
  AC_DEFINE(ThePEG_HAS_HZTOOL, 1, [Set if HZTOOL is present and working])
  AC_SUBST(HZTOOL_LIBS, [], [Libraries needed for HZTOOL])
  AC_SUBST(HZTOOL_LDFLAGS, [], [Link flags needed for HZTOOL])
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi
LIBS="$oldLIB"
LDFLAGS="$oldLDFLAGS"
AM_CONDITIONAL([USEHZTOOL], [test "x$HAS_HZTOOL" == "xyes"])



AC_SUBST([F77],[$FC])
THEPEG_EMPTY_SUBST

AM_CXXFLAGS="-ansi -pedantic -Wall"
AC_SUBST(AM_CXXFLAGS)

AC_CONFIG_FILES([Misc/Makefile
Decay/Makefile
LesHouches/Makefile
HZTOOL/Makefile
Analysis/Makefile
Makefile])
## don't know what to do here  'Doc/Makefile'

AC_OUTPUT

