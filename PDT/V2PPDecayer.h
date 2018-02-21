// -*- C++ -*-
//
// V2PPDecayer.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_V2PPDecayer_H
#define THEPEG_V2PPDecayer_H
// This is the declaration of the V2PPDecayer class.

#include "ThePEG/PDT/FlatDecayer.h"

namespace ThePEG {

/**
 * The V2PPDecayer class performs the decay of a vector meson into two
 * pseudo-scalars according to a flat phase space. If, however the
 * decaying particle comes from a pseudo-scalar and has only one
 * sibling which is a pseudo-scalar (or a photon) the decay is
 * reweighted with \f$\cos^2\f$ (\f$\sin^2\f$ for photon) of the angle
 * between one of the decay products and its grand parent.
 *
 * @see \ref V2PPDecayerInterfaces "The interfaces"
 * defined for V2PPDecayer.
 * @see FlatDecayer
 * @see ParticleData
 */
class V2PPDecayer: public FlatDecayer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Destructor.
   */
  virtual ~V2PPDecayer();
  //@}

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
   * ratios.
   */
  virtual double reweight(const DecayMode & dm, const Particle & parent,
			  const ParticleVector & children) const;
  //@}

public:


  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:


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
   * The grand parent in case reweighting should be done.
   */
  mutable tPPtr grandParent;

  /**
   * The decaying particles sibling in case reweighting should be done.
   */
  mutable tPPtr sibling;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<V2PPDecayer> initV2PPDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  V2PPDecayer & operator=(const V2PPDecayer &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of V2PPDecayer. */
template <>
struct BaseClassTrait<V2PPDecayer,1>: public ClassTraitsType {
  /** Typedef of the first base class of . */
  typedef FlatDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  V2PPDecayer class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<V2PPDecayer>
  : public ClassTraitsBase<V2PPDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::V2PPDecayer"; }
  /** Return the name of the shared library be loaded to get access to
   *  the V2PPDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "V2PPDecayer.so"; }
};

/** @endcond */

}

#endif /* THEPEG_V2PPDecayer_H */
