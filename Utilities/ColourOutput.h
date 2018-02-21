// -*- C++ -*-
//
// ColourOutput.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ColourOutput_H
#define ThePEG_ColourOutput_H

#include <ostream>

namespace ThePEG {

enum class ANSI {
  reset      =  0,

  black      = 30,
  red        = 31,
  green      = 32,
  yellow     = 33,
  blue       = 34,
  magenta    = 35,
  cyan       = 36,
  white      = 37,

  fg_reset   = 39,

  bg_black   = 40,
  bg_red     = 41,
  bg_green   = 42,
  bg_yellow  = 43,
  bg_blue    = 44,
  bg_magenta = 45,
  bg_cyan    = 46,
  bg_white   = 47,

  bg_reset   = 49

};

std::ostream& operator<<(std::ostream & os, ANSI c);

}

#endif /* ThePEG_ColourOutput_H */
