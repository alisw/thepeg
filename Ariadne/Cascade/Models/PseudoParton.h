// -*- C++ -*-
#ifndef Ariadne5_PseudoParton_H
#define Ariadne5_PseudoParton_H
//
// This is the declaration of the PseudoParton class.
//

#include "ThePEG/Config/ThePEG.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * Here is the documentation of the PseudoParton class.
 */
class PseudoParton {

public:

  /**
   * Enum different special partons.
   */
  enum Special { normal = 0,  /**< Normal parton */
		 remnant = 1, /**< Remnant parton */
		 colres = 2   /**< Colour resonance decay product */
  }; 

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  PseudoParton(): isGluon(false), noRecoil(false), type(normal) {}

  /**
   * The standard constructor.
   */
  PseudoParton(bool glue, bool norec, Special intype,
	       const Lorentz5Momentum & pin, tParPtr rp)
    : isGluon(glue), noRecoil(norec), type(intype), p(pin), realParton(rp) {}
  //@}

public:

  /**
   * Is this a quark or a gluon?
   */
  bool isGluon;

  /**
   * Does this parton dislike taking recoil?
   */
  bool noRecoil;

  /**
   * The type of this pseudoparticle.
   */
  Special type;

  /**
   * The momentum of the parton
   */
  Lorentz5Momentum p;

  /**
   * Possibly empty pointer to corresponding real parton.
   */
  tParPtr realParton;

};

}

#endif /* Ariadne5_PseudoParton_H */
