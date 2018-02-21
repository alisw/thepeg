// -*- C++ -*-
//
// HadronizationHandler.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_HadronizationHandler_H
#define ThePEG_HadronizationHandler_H
// This is the declaration of the HadronizationHandler class.

#include "StepHandler.h"

namespace ThePEG {

/**
 * The HadronizationHandler is the base class of all handlers
 * implementing models for hadronization of coloured particles. It is
 * derived from the more general StepHandler class, but does not
 * introduce more functioanality as it stands.
 *
 * @see \ref HadronizationHandlerInterfaces "The interfaces"
 * defined for HadronizationHandler.
 * @see StepHandler
 * @see EventHandler
 * @see SubProcessHandler
 * 
 */
class HadronizationHandler: public StepHandler {

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

private:

  /**
   * Describe an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<HadronizationHandler>
    initHadronizationHandler;

  /**
   *  Private and non-existent assignment operator.
   */
  HadronizationHandler & operator=(const HadronizationHandler &);

};

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of HadronizationHandler.
 */
template <>
struct BaseClassTrait<HadronizationHandler,1>: public ClassTraitsType {
  /** Typedef of the base class of HadronizationHandler. */
  typedef StepHandler NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * HadronizationHandler class.
 */
template <>
struct ClassTraits<HadronizationHandler>:
    public ClassTraitsBase<HadronizationHandler> {
  /** Return the class name. */
  static string className() { return "ThePEG::HadronizationHandler"; }
};

/** @endcond */

}

#endif /* ThePEG_HadronizationHandler_H */
