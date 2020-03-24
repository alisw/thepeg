// -*- C++ -*-
//
// MultipleInteractionHandler.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_MultipleInteractionHandler_H
#define ThePEG_MultipleInteractionHandler_H
// This is the declaration of the MultipleInteractionHandler class.

#include "StepHandler.h"

namespace ThePEG {

/**
 * The MultipleInteractionHandler is the base class of all
 * handlers implementing models for multiple interactions. It is
 * derived from the more general StepHandler class,
 * and does not introduce more functioanality as it stands.
 *
 * @see \ref MultipleInteractionHandlerInterfaces "The interfaces"
 * defined for MultipleInteractionHandler.
 * @see StepHandler
 * @see EventHandler
 * @see SubProcessHandler
 */
class MultipleInteractionHandler: public StepHandler {

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

private:

  /**
   * Describe an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<MultipleInteractionHandler>
    initMultipleInteractionHandler;

  /**
   *  Private and non-existent assignment operator.
   */
  MultipleInteractionHandler & operator=(const MultipleInteractionHandler &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of MultipleInteractionHandler.
 */
template <>
struct BaseClassTrait<MultipleInteractionHandler,1>: public ClassTraitsType {
  /** Typedef of the base class of MultipleInteractionHandler. */
  typedef StepHandler NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * MultipleInteractionHandler class.
 */
template <>
struct ClassTraits<MultipleInteractionHandler>:
    public ClassTraitsBase<MultipleInteractionHandler> {
  /** Return the class name. */
  static string className() { return "ThePEG::MultipleInteractionHandler"; }
};

/** @endcond */

}

#endif /* ThePEG_MultipleInteractionHandler_H */
