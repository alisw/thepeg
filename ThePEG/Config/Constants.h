// -*- C++ -*-
//
// Constants.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Constants_H
#define ThePEG_Constants_H

// This file defines a number of useful constants, placed in the
// namespace <!id>ThePEG::Constants<!!id>.

#include "Unitsystem.h"
#include <cmath>
#include <cfloat>

namespace ThePEG {

/**
 * The Constants namespace containing some useful physical constants
 * with suitable units.
 */
namespace Constants {

using namespace ThePEG::Units;

/** A really large length. */
const Length MaxLength = 1.0e20 * meter;

/** A really large energy. */
const Energy MaxEnergy = 1.0e6 * GeV;

/** A really large squared energy. */
const Energy2 MaxEnergy2 = MaxEnergy * MaxEnergy;

/** The largest possible double. */
const double MaxDouble = DBL_MAX;

/** A really large double. */
const double HugeDouble = DBL_MAX * 1.0e-4;

/** The largest possible float. */
const double MaxFloat = FLT_MAX;

/** A really large floa.t */
const double HugeFloat = FLT_MAX * 0.01;

/** A really large rapidity */
const double MaxRapidity = 100.0;

/** Good old \f$\pi\f$. */
const double pi    = M_PI;

/** Good old \f$2\pi\f$. */
const double twopi = 2.0 * pi;

/** A really large integer */
const long MaxInt = 1000000000L;

/** The smallest non-zero double. */
const double epsilon = DBL_EPSILON;

/** The Euler gamma */
const double EulerGamma = 0.5772156649015329;
}

}

#endif /* ThePEG_Constants_H */
