// -*- C++ -*-
//
// AlphaEMBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_AlphaEMBase_H
#define ThePEG_AlphaEMBase_H
// This is the declaration of the AlphaEMBase class.

#include "RunningCoupling.h"

namespace ThePEG {

/**
 * AlphaEMBase an abstract base class used by the StandardModelBase
 * class to implement the electro-magnetic coupling. Concrete
 * sub-classes must implement the value(Energy2, const
 * StandardModelBase &) function.
 *
 * @see \ref AlphaEMBaseInterfaces "The interfaces"
 * defined for AlphaEMBase.
 * @see StandardModelBase
 */
class AlphaEMBase: public RunningCoupling {

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

private:

  /**
   * Describe an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<AlphaEMBase> initAlphaEMBase;

  /**
   *  Private and non-existent assignment operator.
   */
  AlphaEMBase & operator=(const AlphaEMBase &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of AlphaEMBase. */
template <>
struct BaseClassTrait<AlphaEMBase,1>: public ClassTraitsType {
  /** Typedef of the first base class of AlphaEMBase. */
  typedef RunningCoupling NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  AlphaEMBase class. */
template <>
struct ClassTraits<AlphaEMBase>: public ClassTraitsBase<AlphaEMBase> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::AlphaEMBase"; }
};

/** @endcond */

}

#endif /* ThePEG_AlphaEMBase_H */
