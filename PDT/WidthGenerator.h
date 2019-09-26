// -*- C++ -*-
//
// WidthGenerator.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_WidthGenerator_H
#define ThePEG_WidthGenerator_H
// This is the declaration of the WidthGenerator class.

#include "ThePEG/Config/ThePEG.h"
#include "WidthGenerator.fh"
#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Utilities/Selector.h"
#include "ThePEG/PDT/DecayMode.fh"

namespace ThePEG {

/**
 * WidthGenerator is an abstract base class to be used to encapsulate
 * models for giving the partial decay width of a ParticleData given
 * the specified DecayModes.
 *
 * @see \ref WidthGeneratorInterfaces "The interfaces"
 * defined for WidthGenerator.
 * @see ParticleData
 * @see DecayMode
 */
class WidthGenerator: public Interfaced {

public:

  /** A selector of <code>DecayModes</code>s weighted by their
   *  branching ratios. */
  typedef Selector<tDMPtr> DecayMap;

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return true if this object can be used for the given particle
   * type with the given decay map.
   */
  virtual bool accept(const ParticleData &) const = 0;

  /**
   * Given a Particle, calculate a width.
   */
  Energy width(const Particle &) const;

  /**
   * Given a particle type and a mass of an instance of that particle
   * type, calculate a width.
   */
  virtual Energy width(const ParticleData &, Energy m) const = 0;

  /**
   * Given a particle type and a mass and a width of an instance of
   * that particle type, generate a life time.
   */
  virtual Length lifeTime(const ParticleData &, Energy m, Energy w) const;

  /**
   * Return decay map for the given particle type.
   */
  virtual DecayMap rate(const ParticleData &) const = 0;

  /**
   * Return a decay map for a given Particle instance.
   */
  virtual DecayMap rate(const Particle &);
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
  static AbstractNoPIOClassDescription<WidthGenerator> initWidthGenerator;

  /**
   *  Private and non-existent assignment operator.
   */
  WidthGenerator & operator=(const WidthGenerator &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of WidthGenerator. */
template <>
struct BaseClassTrait<WidthGenerator,1>: public ClassTraitsType {
  /** Typedef of the first base class of WidthGenerator. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  WidthGenerator class. */
template <>
struct ClassTraits<WidthGenerator>: public ClassTraitsBase<WidthGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::WidthGenerator"; }
};

/** @endcond */

}

#endif /* ThePEG_WidthGenerator_H */
