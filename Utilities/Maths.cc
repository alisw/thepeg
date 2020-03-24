// -*- C++ -*-
//
// Maths.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include <config.h>
#include "Maths.h"
#include <cmath>

double ThePEG::Math::exp1m(double x) {
  using namespace std;
#ifndef ThePEG_HAS_EXPM1
  return 2.0*std::exp(x/2.0)*std::sinh(-x/2.0);
#else
  return -expm1(x);
#endif
}

double ThePEG::Math::log1m(double x) {
  using namespace std;
#ifndef ThePEG_HAS_LOG1P
#ifndef ThePEG_HAS_ATANH
  volatile double y = 1.0 - x;
  return log(y) - ((y - 1.0) + x)/y ; /* cancels errors with IEEE arithmetic */
#else
  return 2.0*atanh(x/(x-2.0));
#endif
#else
  return log1p(-x);
#endif
}

double ThePEG::Math::powi(double x, int i) {
  switch ( i ) {
  case 0: return 1.0;
  case -1: return 1/x;
  case -2: return 1/x/x;
  case -3: return 1/x/x/x;
  case 1: return x;
  case 2: return x*x;
  case 3: return x*x*x;
  default: return i > 0? powi(x, i - 1)*x: powi(x, i + 1)/x;
  }
}
