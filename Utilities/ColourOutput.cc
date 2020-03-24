// -*- C++ -*-
//
// ColourOutput.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourOutput class.
//

#include "ColourOutput.h"
#include <iostream>
#include <unistd.h>


namespace {

bool useColours(const std::ostream & os) {

	// are we cout?
	if ( os.rdbuf() == std::cout.rdbuf() ) {
		return isatty(fileno(stdout));
	}
	// are we cerr?
	else if ( os.rdbuf() == std::cerr.rdbuf() ) {
		return isatty(fileno(stderr));
	}

	// otherwise play safe
	return false;

}

}

namespace ThePEG {

std::ostream& operator<<(std::ostream & os, ANSI c) {
	if ( useColours(os) )
    os << "\x1B[" << static_cast<int>(c) << 'm';
  return os;
}

}
