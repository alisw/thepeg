// -*- C++ -*-
//
// HelicityDefinitions.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_HelicityDefinitions_H
#define THEPEG_HelicityDefinitions_H
// This is the declaration of the HelicityDefinitions class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Utilities/Exception.h"

/** \file HelicityDefinitions.h
 *
 * This file contains enumerations used by LorentzSpinor and
 * LorentzSpinorBar classes.
 *
 * @see LorentzSpinor
 *
 * @author Peter Richardson
 */

namespace ThePEG {

/**
 * The Helicity namespace contains classes for spin representation
 * classes in ThePEG.
 */
namespace Helicity {

/**
 * Enumeration to specify spinor type.
 */
enum class SpinorType {
  u, /**< u spinor. */
  v, /**< v spinor. */
  unknown /**< Undefined spinor type. */
};

/** \ingroup Helicity
 * Definition of the enumerated values of the phase to include in the 
 * calculation of the polarization vector.
 */
enum VectorPhase {
  vector_phase, /**< Include the phase factor.*/
  vector_nophase, /**< No phase-factor. */
  default_vector_phase=vector_nophase /**< Default option.*/
};
  
/** @cond EXCEPTIONCLASSES */
/** Exception class used by Helicity classes to signal a logical error. */
class HelicityLogicalError: public Exception {};

/** Exception class used by Helicity classes to signal a inconsistencies. */
class HelicityConsistencyError: public Exception {};
/** @endcond */

}
}

#endif /* THEPEG_HelicityDefinitions_H */
