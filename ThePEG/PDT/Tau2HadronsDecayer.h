// -*- C++ -*-
//
// Tau2HadronsDecayer.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_Tau2HadronsDecayer_H
#define THEPEG_Tau2HadronsDecayer_H
// This is the declaration of the Tau2HadronsDecayer class.

#include "ThePEG/PDT/FlatDecayer.h"

namespace ThePEG {

/**
 * The Tau2HadronsDecayer class inherits FlatDecayer and can perform
 * the decays of tau to neutrino + hadrons according to phase space,
 * with an extra weight \f$x_\nu(3-x_\nu)\f$.
 *
 * @see \ref Tau2HadronsDecayerInterfaces "The interfaces"
 * defined for Tau2HadronsDecayer.
 */
class Tau2HadronsDecayer: public FlatDecayer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Destructor.
   */
  virtual ~Tau2HadronsDecayer();
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
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<Tau2HadronsDecayer> initTau2HadronsDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  Tau2HadronsDecayer & operator=(const Tau2HadronsDecayer &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of Tau2HadronsDecayer. */
template <>
struct BaseClassTrait<Tau2HadronsDecayer,1>: public ClassTraitsType {
  /** Typedef of the first base class of Tau2HadronsDecayer. */
  typedef FlatDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  Tau2HadronsDecayer class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<Tau2HadronsDecayer>
  : public ClassTraitsBase<Tau2HadronsDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::Tau2HadronsDecayer"; }
  /** Return the name of the shared library be loaded to get access to
   *  the Tau2HadronsDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "Tau2HadronsDecayer.so"; }
};

/** @endcond */

}

#endif /* THEPEG_Tau2HadronsDecayer_H */
