// -*- C++ -*-
//
// Onium3GDecayer.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_Onium3GDecayer_H
#define THEPEG_Onium3GDecayer_H
// This is the declaration of the Onium3GDecayer class.

#include "ThePEG/PDT/FlatDecayer.h"

namespace ThePEG {

/**
 * The Onium3GDecayer class inherits from performs FlatDecayer and
 * will reweight the flat phase space suitable to describe the decay
 * of a spin-1 onium resonance into three gluons or two gluons and a
 * photon. After the decay the collision handler is instructed to
 * restart the generation from the hadronization (or optionally the
 * parton cascade) stage.
 *
 * @see \ref Onium3GDecayerInterfaces "The interfaces"
 * defined for Onium3GDecayer.
 * @see FlatDecayer
 * @see ParticleData
 */
class Onium3GDecayer: public FlatDecayer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  Onium3GDecayer() : doShower(true), theMinGGMass(2.0*GeV) {}

  /**
   * Destructor.
   */
  virtual ~Onium3GDecayer();
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

  /**
   * Return true if the produced gluons should be showered.
   */
  bool shower() const { return doShower; }

  /**
   * Return the minimum invariant mass between two gluons in gamma-g-g
   * decays.
   */
  Energy minGGMass() const { return theMinGGMass; }

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
   * If true the produced gluons should be showered.
   */
  bool doShower;

  /**
   * The minimum invariant mass between two gluons in gamma-g-g
   * decays.
   */
  Energy theMinGGMass;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<Onium3GDecayer> initOnium3GDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  Onium3GDecayer & operator=(const Onium3GDecayer &) = delete;

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of Onium3GDecayer. */
template <>
struct BaseClassTrait<Onium3GDecayer,1>: public ClassTraitsType {
  /** Typedef of the first base class of Onium3GDecayer. */
  typedef FlatDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  Onium3GDecayer class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<Onium3GDecayer>
  : public ClassTraitsBase<Onium3GDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::Onium3GDecayer"; }
  /** Return the name of the shared library be loaded to get access to
   *  the Onium3GDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "Onium3GDecayer.so"; }
};

/** @endcond */

}

#endif /* THEPEG_Onium3GDecayer_H */
