// -*- C++ -*-
//
// FlatDecayer.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_FlatDecayer_H
#define ThePEG_FlatDecayer_H
// This is the declaration of the FlatDecayer class.

#include "ThePEG/PDT/Decayer.h"

namespace ThePEG {

/**
 * The FlatDecayer class inrerits from the abstract Decayer class and
 * implements the decay of a given Particle to a given set of children
 * according to a flat phase space distribution.
 *
 * It is possible to implement a more complicated decay distribution
 * by inheriting from the FlatDecayer class and only override the
 * virtual function reweight() to return a weight (between zero and
 * one) of a given phase space point relative to the flat
 * distribution.
 *
 * @see \ref FlatDecayerInterfaces "The interfaces"
 * defined for FlatDecayer.
 * @see Decayer
 * @see Particle
 */
class FlatDecayer: public Decayer {

public:

  /** @name Virtual functions required by the Decayer class.
   */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;

  /**
   * Give a weight to a phase space point. To be overridden by
   * subclasses. For a given decay mode, \a dm, decaying \a parent
   * particle and decayproducts, \a children, distributed according to
   * a flat distribution in phase space, return a weight (less or
   * equal to unity) modifying the flat distribution to the desired
   * one. Note that the chosen phase space point may be rejected, but
   * the chosen decay channel will not. This means that the weight
   * returned by this function does not influence the branching
   * ratios. For the FlatDecayer class this function simply returns 1.
   */
  virtual double reweight(const DecayMode &, const Particle & ,
			  const ParticleVector & ) const {
    return 1.0;
  }
  //@}

public:

  /**
   * Standard Init function used to initialize the interfaces.
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
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<FlatDecayer> initFlatDecayer;

  /**
   *  Private and non-existent assignment operator.
   */
  FlatDecayer & operator=(const FlatDecayer &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of FlatDecayer. */
template <>
struct BaseClassTrait<FlatDecayer,1>: public ClassTraitsType {
  /** Typedef of the first base class of FlatDecayer. */
  typedef Decayer NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  FlatDecayer class. */
template <>
struct ClassTraits<FlatDecayer>: public ClassTraitsBase<FlatDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::FlatDecayer"; }
};

/** @endcond */

}

#endif /* ThePEG_FlatDecayer_H */
