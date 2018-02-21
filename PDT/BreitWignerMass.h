// -*- C++ -*-
//
// BreitWignerMass.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_BreitWignerMass_H
#define ThePEG_BreitWignerMass_H
// This is the declaration of the BreitWignerMass class.

#include "ThePEG/PDT/MassGenerator.h"

namespace ThePEG {

/**
 * BreitWignerMass is derived from MassGenerator and is able to
 * generate the mass for a particle given its nominal mass and its
 * with.
 *
 *
 * @see MassGenerator
 * @see ParticleData
 * 
 */
class BreitWignerMass: public MassGenerator {

public:

  /** @name Virtual methods required by the MassGenerator base class. */
  //@{
  /**
   * Return true if this mass generator can handle the given particle
   * type.
   */
  virtual bool accept(const ParticleData &) const { return true; }

  /**
   * Generate a mass for an instance of a given particle type.
   */
  virtual Energy mass(const ParticleData &) const;
  //@}

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * Describe concrete class without persistent data.
   */
  static NoPIOClassDescription<BreitWignerMass> initBreitWignerMass;

  /**
   *  Private and non-existent assignment operator.
   */
  BreitWignerMass & operator=(const BreitWignerMass &);

};


/** @cond TRAITSPECIALIZATIONS */
ThePEG_DECLARE_DYNAMIC_CLASS_TRAITS(BreitWignerMass,MassGenerator,"BreitWignerMass.so");
/** @endcond */

}

#endif /* ThePEG_BreitWignerMass_H */
