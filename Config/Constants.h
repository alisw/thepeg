// -*- C++ -*-
//
// Constants.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
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
constexpr Length MaxLength = 1.0e23_mm;

/** A really large energy. */
constexpr Energy MaxEnergy = 1.0e6_GeV;

/** A really large squared energy. */
constexpr Energy2 MaxEnergy2 = MaxEnergy * MaxEnergy;

/** The largest possible double. */
constexpr double MaxDouble = DBL_MAX;

/** A really large double. */
constexpr double HugeDouble = DBL_MAX * 1.0e-4;

/** The largest possible float. */
constexpr double MaxFloat = FLT_MAX;

/** A really large floa.t */
constexpr double HugeFloat = FLT_MAX * 0.01;

/** A really large rapidity */
constexpr double MaxRapidity = 100.0;

/** Good old \f$\pi\f$. */
constexpr double pi    = M_PI;

/** Good old \f$2\pi\f$. */
constexpr double twopi = 2.0 * pi;

/** A really large integer */
constexpr long MaxInt = 1000000000L;

/** The smallest non-zero double. */
constexpr double epsilon = DBL_EPSILON;

/** The Euler gamma */
constexpr double EulerGamma = 0.5772156649015329;

/** \f$\zeta(2)\f$. */
constexpr double zeta2 = pi*pi/6.;

/** \f$\zeta(3)\f$. */
constexpr double zeta3 = 1.2020569031595943;

/** \f$\zeta(4)\f$. */
constexpr double zeta4 = 0.4*zeta2*zeta2;

/** \f$\zeta(5)\f$. */
constexpr double zeta5 = 1.0369277551433699;

/** \f$\zeta(6)\f$. */
constexpr double zeta6 = 4.*zeta2*zeta4/7.;
}

}

#endif /* ThePEG_Constants_H */
