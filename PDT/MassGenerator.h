// -*- C++ -*-
//
// MassGenerator.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_MassGenerator_H
#define ThePEG_MassGenerator_H
// This is the declaration of the MassGenerator class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Interface/Interfaced.h"
#include "MassGenerator.fh"

namespace ThePEG {

/**
 * MassGenerator is an abstract base class used by the ParticleData
 * class to generate a mass for a Particle instance.
 *
 * @see \ref MassGeneratorInterfaces "The interfaces"
 * defined for MassGenerator.
 * @see ParticleData
 */
class MassGenerator: public Interfaced {

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return true if this mass generator can handle the given particle
   * type.
   */
  virtual bool accept(const ParticleData &) const = 0;

  /**
   * Generate a mass for an instance of a given particle type.
   */
  virtual Energy mass(const ParticleData &) const = 0;
  //@}

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<MassGenerator> initMassGenerator;

  /**
   *  Private and non-existent assignment operator.
   */
  MassGenerator & operator=(const MassGenerator &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of MassGenerator. */
template <>
struct BaseClassTrait<MassGenerator,1>: public ClassTraitsType {
  /** Typedef of the first base class of MassGenerator. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  MassGenerator class. */
template <>
struct ClassTraits<MassGenerator>: public ClassTraitsBase<MassGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::MassGenerator"; }
};

/** @endcond */

}

#endif /* ThePEG_MassGenerator_H */
