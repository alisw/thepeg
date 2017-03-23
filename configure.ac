dnl Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([ThePEG],[2.0.4],[http://www.thep.lu.se/ThePEG/],[ThePEG])
AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_MACRO_DIR([m4])

THEPEG_LIBTOOL_VERSION_INFO(20,0,0)

AC_CONFIG_SRCDIR([EventRecord/SubProcess.h])
AC_CONFIG_HEADERS([Config/config.h Config/ThePEG_Qty.h:Config/ThePEG_Qty.h.a:Config/ThePEG_Qty.h.in:Config/ThePEG_Qty.h.b])

AC_CANONICAL_HOST

case "${host}" in
  *-darwin[[0156]].*)
    AC_MSG_ERROR([ThePEG requires OS X 10.3 or later])
    ;;
  *-darwin7.*)
    if test "x$MACOSX_DEPLOYMENT_TARGET" != "x10.3"; then
       AC_MSG_ERROR(
       [Please export the MACOSX_DEPLOYMENT_TARGET variable, set to 10.3])
    fi
    ;;
esac

AC_LANG(C++)
AM_INIT_AUTOMAKE([1.9 gnu dist-bzip2 subdir-objects -Wall])
dnl also include std-options once --version and --help exist
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

dnl Checks for C++ compiler. Handle C++11 flags.
AC_PROG_CXX

dnl this block can go once we decide to have C++11 always on
AC_MSG_CHECKING([whether to include C++11 flag for testing])
AC_ARG_ENABLE(stdcxx11,
        AC_HELP_STRING([--enable-stdcxx11],
    [turn on C++11 flag (only for testing, do not use in production!)]),
        [],
        [enable_stdcxx11=no]
        )


if test "x$enable_stdcxx11" = "xyes"; then
AC_MSG_RESULT([yes])
dnl remove the wrapper if block and "optional" once we decide to have C++11 always on
AX_CXX_COMPILE_STDCXX([11],[noext],[optional])
if test "x$HAVE_CXX11" != "x1"; then
AC_MSG_ERROR([compiler does not recognize requested c++11 option])
fi
else
AC_MSG_RESULT([no])
fi

dnl check for unistd
AC_CHECK_HEADER([unistd.h],[],
      [AC_MSG_ERROR([ThePEG needs "unistd.h". Check your system is POSIX-compliant.])])

dnl Checks for programs.
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_LN_S

LT_PREREQ([2.2])
LT_INIT([disable-static dlopen pic-only])

VL_LIB_READLINE

THEPEG_UNIT_CHECKING
THEPEG_CHECK_GSL
THEPEG_SEARCH_LHAPDF
THEPEG_CHECK_HEPMC
THEPEG_CHECK_RIVET
FASTJET_CHECK_FASTJET

AX_CHECK_ZLIB

dnl AX_CHECK_BZ2LIB

THEPEG_DEFINE_ENVDEFAULT(ThePEG_GZREAD_FILE,GZREAD_FILE,gunzip -c,[The command which, taking the name of a gzipped file as argument, unzips it and prints it to stdout. Default is "gunzip -c"])

THEPEG_DEFINE_ENVDEFAULT(ThePEG_GZWRITE_FILE,GZWRITE_FILE,[gzip -c > ],[The command which, taking the name of a gzipped file as argument, reads stdin, zips it and writes it to the file. Default is "gzip -c > ".])

THEPEG_DEFINE_ENVDEFAULT(ThePEG_BZ2READ_FILE,BZ2READ_FILE,bunzip2 -c,[The command which, taking the name of a bzipped file as argument, unzips it and prints it to stdout. Default is "bunzip2 -c".])

THEPEG_DEFINE_ENVDEFAULT(ThePEG_BZ2WRITE_FILE,BZ2WRITE_FILE,[bzip2 -c > ],[The command which, taking the name of a bzipped file as argument, reads stdin, zips it and writes it to the file. Default is "bzip2 -c > ".])

THEPEG_CHECK_EXPM1
THEPEG_CHECK_LOG1P
THEPEG_CHECK_DLOPEN

AX_COMPILER_VENDOR
case "${ax_cv_cxx_compiler_vendor}" in
    clang)
        case "${host}" in
            *linux*)
                AC_MSG_WARN([
*****************************************************************************
 clang/LLVM ignores the CPU floating-point environment.
 All floating point exception trapping will be disabled.
*****************************************************************************])
                ;;
        esac
esac
THEPEG_CHECK_FPUCONTROL
THEPEG_CHECK_FENV

AM_CPPFLAGS="-I\$(top_builddir)/include \$(GSLINCLUDE)"

case "${ax_cv_cxx_compiler_vendor}" in
     gnu)
    AM_CXXFLAGS="-pedantic -Wall -W"
    ;;
     clang)
    AM_CXXFLAGS="-pedantic -Wall -Wno-overloaded-virtual -Wno-unused-function"
dnl  -Wno-unneeded-internal-declaration
    ;;
     intel)
    AM_CXXFLAGS="-strict-ansi -Wall -wd13000,1418,981,444,383,1599,1572,2259,980"
    ;;
esac

THEPEG_CHECK_ABS_BUG

AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)

dnl do an actual capability check on ld instead of this workaround
case "${host}" in
  *-darwin*)
     ;;
  *)
     AM_LDFLAGS="-Wl,--enable-new-dtags"
     ;;
esac

AC_SUBST(AM_LDFLAGS)

THEPEG_EMPTY_SUBST

AC_PATH_PROG(PERL, perl)

AC_ARG_WITH(javagui,
            [  --with-javagui          Compile and install the java-based GUI.])

if test "x$with_javagui" != "xno"; then
  THEPEG_HAS_JAVA([1.4], [], [with_javagui=no; AC_MSG_NOTICE([Java GUI disabled])])
fi

AM_CONDITIONAL([JAVAGUI], [test "x$with_javagui" != "xno"])

AC_CONFIG_FILES([Helicity/Makefile
Helicity/WaveFunction/Makefile
Helicity/Vertex/Makefile
Helicity/Vertex/Scalar/Makefile
Helicity/Vertex/Vector/Makefile
Helicity/Vertex/Tensor/Makefile
Utilities/Makefile
include/Makefile
Interface/Makefile
LesHouches/Makefile
Vectors/Makefile
PDT/Makefile
PDF/Makefile
Persistency/Makefile
Config/Makefile
Handlers/Makefile
MatrixElement/Makefile
Pointer/Makefile
lib/Makefile
lib/Makefile.common.install
src/Makefile
ACDC/Makefile
Repository/Makefile
EventRecord/Makefile
StandardModel/Makefile
Cuts/Makefile
Analysis/Makefile
Doc/Makefile
Doc/MakeDocs.in
Doc/refman.h
Doc/refman.conf
java/Makefile
Makefile])
AC_CONFIG_FILES([Doc/fixinterfaces.pl],[chmod +x Doc/fixinterfaces.pl])

BOOST_REQUIRE([1.41],[AC_SUBST([BOOST_NOT_FOUND],[true])])
BOOST_TEST()
THEPEG_BOOST_UNIT_TEST

THEPEG_OVERVIEW

AC_CONFIG_COMMANDS([summary],[cat config.thepeg])

AC_OUTPUT
