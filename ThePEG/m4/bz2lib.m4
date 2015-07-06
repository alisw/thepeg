# ===========================================================================
#          http://www.nongnu.org/autoconf-archive/ax_check_zlib.html
# modified for bz2lib 
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_BZ2LIB()
#
# DESCRIPTION
#
#   This macro searches for an installed bz2lib library. If nothing was
#   specified when calling configure, it searches first in /usr/local and
#   then in /usr. If the --with-bz2lib=DIR is specified, it will try to find
#   it in DIR/include/bzlib.h and DIR/lib/libbz2.a. If --without-bz2lib is
#   specified, the library is not searched at all.
#
#   If either the header file (bzlib.h) or the library (libbz2) is not found,
#   the configuration exits on error, asking for a valid bz2lib installation
#   directory or --without-bz2lib.
#
#   The macro defines the symbol HAVE_LIBBZ2 if the library is found. You
#   should use autoheader to include a definition for this symbol in a
#   config.h file. Sample usage in a C/C++ source is as follows:
#
#     #ifdef HAVE_LIBBZ2
#     #include <bzlib.h>
#     #endif /* HAVE_LIBBZ2 */
#
# LICENSE
#
#   Copyright (c) 2008 Loic Dachary <loic@senga.org>
#   changes for BZ2 and ThePEG: Copyright (c) 2010 David Grellscheid
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 4

AU_ALIAS([CHECK_BZ2LIB], [AX_CHECK_BZ2LIB])
AC_DEFUN([AX_CHECK_BZ2LIB],
#
# Handle user hints
#
[AC_MSG_CHECKING(if bz2lib is wanted)
AC_ARG_WITH(bz2lib,
[  --with-bz2lib=DIR root directory path of bz2lib installation [defaults to
                    /usr/local or /usr if not found in /usr/local]
  --without-bz2lib to disable bz2lib usage completely],
[if test "$withval" != no ; then
  AC_MSG_RESULT(yes)
  if test -d "$withval"
  then
    BZ2LIB_HOME="$withval"
  else
    AC_MSG_WARN([Sorry, $withval does not exist, checking usual places])
  fi
else
  AC_MSG_RESULT(no)
fi])

BZ2LIB_HOME=/usr/local
if test ! -f "${BZ2LIB_HOME}/include/bzlib.h"
then
        BZ2LIB_HOME=/usr
fi

#
# Locate bz2lib, if wanted
#
if test -n "${BZ2LIB_HOME}"
then
        BZ2LIB_OLD_LDFLAGS=$LDFLAGS
        BZ2LIB_OLD_CPPFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS -L${BZ2LIB_HOME}/lib"
        CPPFLAGS="$CPPFLAGS -I${BZ2LIB_HOME}/include"
        AC_LANG_SAVE
        AC_LANG_C
        AC_CHECK_LIB(bz2, BZ2_bzDecompressEnd, [bz2lib_cv_libbz2=yes], [bz2lib_cv_libbz2=no])
        AC_CHECK_HEADER(bzlib.h, [bz2lib_cv_bz2lib_h=yes], [bz2lib_cv_bz2lib_h=no])
        AC_LANG_RESTORE
        if test "$bz2lib_cv_libbz2" = "yes" -a "$bz2lib_cv_bz2lib_h" = "yes"
        then
                #
                # If both library and header were found, use them
                #
                AC_CHECK_LIB(bz2, BZ2_bzDecompressEnd)
                AC_MSG_CHECKING(bz2lib in ${BZ2LIB_HOME})
                AC_MSG_RESULT(ok)
        else
                #
                # If either header or library was not found, revert and bomb
                #
                AC_MSG_CHECKING(bz2lib in ${BZ2LIB_HOME})
                LDFLAGS="$BZ2LIB_OLD_LDFLAGS"
                CPPFLAGS="$BZ2LIB_OLD_CPPFLAGS"
                AC_MSG_RESULT(failed)
                AC_MSG_ERROR(either specify a valid bz2lib installation with --with-bz2lib=DIR or disable bz2lib usage with --without-bz2lib)
        fi
fi

])
