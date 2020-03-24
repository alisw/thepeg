// -*- C++ -*-
//
// PtGenerator.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_PtGenerator_H
#define ThePEG_PtGenerator_H
// This is the declaration of the PtGenerator class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Handlers/HandlerBase.h"

namespace ThePEG {

/**
 * PtGenerator is the base for all classes implementing alternative
 * models for transverse momentum generation.  It inherits from the
 * HandlerBase which among other things provides forward access to the
 * random number object held by the EventGenerator object.
 *
 * @see \ref PtGeneratorInterfaces "The interfaces"
 * defined for PtGenerator.
 * @see HandlerBase
 * @see EventGenerator
 */
class PtGenerator: public HandlerBase {

public:

  /** @name Virtual functions to be implemented by sub-classes. */
  //@{
  /**
   * Generate (\f$k_x, k_y\f$) components of the transverse
   * momentum.
   */
  virtual TransverseMomentum generate() const =0;
  //@}

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

private:

  /**
   * Describe an abstract class without persistent data.
   */
  static AbstractClassDescription<PtGenerator> initPtGenerator;

  /**
   * Private and non-existent assignment operator.
   */
   PtGenerator & operator=(const PtGenerator &) = delete;

};


/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of PtGenerator.
 */
template <>
struct BaseClassTrait<PtGenerator,1>: public ClassTraitsType {
  /** Typedef of the base class of PtGenerator. */
  typedef HandlerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * PtGenerator class.
 */
template <>
struct ClassTraits<PtGenerator>: public ClassTraitsBase<PtGenerator> {
  /** Return the class name. */
  static string className() { return "ThePEG::PtGenerator"; }
};


/** @endcond */

}

#endif /* ThePEG_PtGenerator_H */
