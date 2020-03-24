// -*- C++ -*-
//
// Complex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Complex_H
#define ThePEG_Complex_H
//
// This is file wraps the standard complex header and makes some
// convenient typedefs in the ThePEG namespace.
//

#include <complex>

namespace ThePEG {

using std::complex;

/// ThePEG code should use Complex for all complex scalars
using Complex = std::complex<double>;

}

#endif /* ThePEG_Complex_H */

