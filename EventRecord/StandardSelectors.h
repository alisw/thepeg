// -*- C++ -*-
//
// StandardSelectors.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_StandardSelectors_H
#define ThePEG_StandardSelectors_H

/**
 * @file StandardSelectors.h This file contains declarations of
 * standard selector classes. The classes contain only static
 * functions and are assumed to be used as template arguments to the
 * ParticleSelector class.
 */

#include "SelectorBase.h"
#include "ParticleTraits.h"

namespace ThePEG {

/**
 * The AllSelector class is used to extract all particles from an event.
 *
 * @see SelectorBase
 * @see ParticleSelector
 */
struct AllSelector: public SelectorBase {

  /**
   * Static method corresponding to the virtual check()
   * method. Returns true.
   */
  static bool Check(const Particle &) { return true; }

  /**
   * Static method corresponding to the virtual intermediate()
   * method. Returns true.
   */
  static bool Intermediate() { return true; }

  /**
   * Static method corresponding to the virtual finalState()
   * method. Returns true.
   */
  static bool FinalState() { return true; }

  /**
   * Static method corresponding to the virtual allSteps()
   * method. Returns true.
   */
  static bool AllSteps() { return true; }

  /**
   * Static method corresponding to the virtual allCollisions()
   * method. Returns true.
   */
  static bool AllCollisions() { return true; }

};

/** Typedef to declare a selector used to extract all particles from an
 *  event. */
typedef ParticleSelector<AllSelector> SelectAll;


/**
 * The FinalStateSelector class is used to extract all final state particles
 * from an event.
 *
 * @see SelectorBase
 * @see ParticleSelector
 */
struct FinalStateSelector: public SelectorBase {

  /**
   * Static method corresponding to the virtual intermediate()
   * method. Returns false.
   */
  static bool Intermediate() { return false; }

  /**
   * Static method corresponding to the virtual allSteps()
   * method. Returns false.
   */
  static bool AllSteps() { return false; }

};

/** Typedef to declare a selector used to extract all final state
 *  particles from an event. */
typedef ParticleSelector<FinalStateSelector> SelectFinalState;

/**
 * The IntermediateSelector class is used to extract only intermediate
 * particles from an event.
 *
 * @see SelectorBase
 * @see ParticleSelector
 */
struct IntermediateSelector: public SelectorBase {

  /**
   * Static method corresponding to the virtual check()
   * method. Returns true.
   */
  static bool Check(const Particle &) { return true; }

  /**
   * Static method corresponding to the virtual intermediate()
   * method. Returns true.
   */
   static bool Intermediate() { return true; }

   /**
   * Static method corresponding to the virtual finalState()
   * method. Returns false.
   */
  static bool FinalState() { return false; }

  /**
   * Static method corresponding to the virtual allSteps()
   * method. Returns true.
   */
  static bool AllSteps() { return true; }

  /**
   * Static method corresponding to the virtual allCollisions()
   * method. Returns true.
   */
  static bool AllCollisions() { return true; }

};

/** Typedef to declare a selector used to extract all intermediate
 *  particles from an event. */
typedef ParticleSelector<IntermediateSelector> SelectIntermediates;

/**
 * The PrimaryCollisionSelector class is used to extract all particles
 * from the primary Collision of an event.
 *
 * @see SelectorBase
 * @see ParticleSelector
 */
struct PrimaryCollisionSelector: public SelectorBase {

  /**
   * Static method corresponding to the virtual check()
   * method. Returns true.
   */
  static bool Check(const Particle &) { return true; }

  /**
   * Static method corresponding to the virtual intermediate()
   * method. Returns true.
   */
  static bool Intermediate() { return true; }

  /**
   * Static method corresponding to the virtual finalState()
   * method. Returns true.
   */
  static bool FinalState() { return true; }

  /**
   * Static method corresponding to the virtual allSteps()
   * method. Returns true.
   */
  static bool AllSteps() { return true; }

  /**
   * Static method corresponding to the virtual allCollisions()
   * method. Returns false.
   */
  static bool AllCollisions() { return false; }

};

/** Typedef to declare a selector used to extract all particles from
 *  the primary Collision of an event. */
typedef ParticleSelector<PrimaryCollisionSelector> SelectPrimaryCollision;

/**
 * The ChargedSelector class is used to extract all charged particles
 * from an event.
 *
 * @see SelectorBase
 * @see ParticleSelector
 */
struct ChargedSelector: public SelectorBase {

  /**
   * Static method corresponding to the virtual check()
   * method. Returns true if the given particle is charged..
   */
  static bool Check(const Particle & p) {
    return ParticleTraits<Particle>::iCharge(p);
  }

};

/** Typedef to declare a selector used to extract all charged
 *  particles from an event. */
typedef ParticleSelector<ChargedSelector> SelectCharged;

}

#endif /* ThePEG_StandardSelectors_H */
