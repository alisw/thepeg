dnl http://autoconf-archive.cryp.to/vl_lib_readline.html
dnl
dnl Copyright © 2008 Ville Laurikari <vl@iki.fi>
dnl
dnl Modifications for ThePEG: 
dnl Copyright © 2009 D.Grellscheid <herwig@projects.hepforge.org>
dnl
dnl Copying and distribution of this file, with or without modification, 
dnl are permitted in any medium without royalty provided the copyright 
dnl notice and this notice are preserved. 
dnl

AC_DEFUN([VL_LIB_READLINE], [
  AC_ARG_ENABLE(readline,
  AC_HELP_STRING([--disable-readline],[turns off readline support.]),
        [],
        [enable_readline=yes]
        )

  if test "$enable_readline" = "yes"; then

  AC_CACHE_CHECK([for a readline compatible library],
                 vl_cv_lib_readline, [
    ORIG_LIBS="$LIBS"
    for readline_lib in readline edit editline; do
      for termcap_lib in "" termcap curses ncurses; do
        if test -z "$termcap_lib"; then
          TRY_LIB="-l$readline_lib"
        else
          TRY_LIB="-l$readline_lib -l$termcap_lib"
        fi
        LIBS="$ORIG_LIBS $TRY_LIB"
        AC_TRY_LINK_FUNC(readline, vl_cv_lib_readline="$TRY_LIB")
        if test -n "$vl_cv_lib_readline"; then
          break
        fi
      done
      if test -n "$vl_cv_lib_readline"; then
        break
      fi
    done
    if test -z "$vl_cv_lib_readline"; then
      vl_cv_lib_readline="no"
      LIBS="$ORIG_LIBS"
    fi
  ])

  if test "$vl_cv_lib_readline" != "no"; then
    AC_DEFINE(HAVE_LIBREADLINE, 1,
              [Define if you have a readline compatible library])
    AC_CHECK_HEADERS(readline.h readline/readline.h)
    AC_CACHE_CHECK([whether readline supports history],
                   vl_cv_lib_readline_history, [
      vl_cv_lib_readline_history="no"
      AC_TRY_LINK_FUNC(add_history, vl_cv_lib_readline_history="yes")
    ])
    if test "$vl_cv_lib_readline_history" = "yes"; then
      AC_DEFINE(HAVE_READLINE_HISTORY, 1,
                [Define if your readline library has \`add_history'])
      AC_CHECK_HEADERS(history.h readline/history.h)
    fi
  fi

  fi
])dnl
