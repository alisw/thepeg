// -*- C++ -*-
//
// ACDCGenConfig.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ACDCGenConfig_H
#define ACDCGenConfig_H

/** @file ACDCGenConfig.h is the main config header file for
 *  ACDCGen. Do not make changes in this file. If you need to modify
 *  anything, edit a copy of the file which can be included instead of
 *  this file using the macro <code>ACDC_ALT_CONFIG</code>.
 * 
 * ACDCGen uses some classes and functions from the standard
 * library. These are here imported into the ACDCGenerator
 * namespace. If alternative classes with the same API are needed
 * these should be imported with the same name into the namespace in
 * the <code>ACDC_ALT_CONFIG</code> file.
*/

#ifndef ACDC_ALT_CONFIG

#include <vector>
#include <map>
#include <limits>

/** The namespace in which all ACDCGen classes are defined. */
namespace ACDCGenerator {


using std::vector;
using std::multimap;
using std::numeric_limits;
using std::map;
using std::max;
using std::min;
using std::swap;
using std::make_pair;

/** The integer type used to represent the dimension of the a
 *  functions to be sampled- */
typedef short DimType;

}

#else

#include ACDC_ALT_CONFIG

#endif

namespace ACDCGenerator {

/** A vector of doubles. */
typedef vector<double> DVector;

}

#endif
