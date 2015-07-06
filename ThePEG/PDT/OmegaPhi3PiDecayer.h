// -*- C++ -*-
//
// OmegaPhi3PiDecayer.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_OmegaPhi3PiDecayer_H
#define THEPEG_OmegaPhi3PiDecayer_H
// This is the declaration of the OmegaPhi3PiDecayer class.

#include "ThePEG/PDT/FlatDecayer.h"

namespace ThePEG {

/**
 * The OmegaPhi3PiDecayer class inherits from performs FlatDecayer and
 * will reweight the flat phase space suitable to describe the decay
 * of a \f$\phi\f$ or an \f$\omega\f$ into \f$\pi^+\pi^-\pi^0\f$. It
 * will in fact decay anything into \f$\pi^+\pi^-\pi^0\f$ assuming the
 * same matrix element.
 *
 * @see \ref OmegaPhi3PiDecayerInterfaces "The interfaces"
 * defined for OmegaPhi3PiDecayer.
 * @see FlatDecayer
 * @see ParticleData
 */
class OmegaPhi3PiDecayer: public FlatDecayer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  OmegaPhi3PiDecayer() :  margin(150.0) {}
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
   * Used to multiply the bare weight to get something below unity. In
   * the Fortran pythia version it was set to 150 for unknown reasons.
   */
  double margin;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<OmegaPhi3PiDecayer> initOmegaPhi3PiDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  OmegaPhi3PiDecayer & operator=(const OmegaPhi3PiDecayer &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of OmegaPhi3PiDecayer. */
template <>
struct BaseClassTrait<OmegaPhi3PiDecayer,1>: public ClassTraitsType {
  /** Typedef of the first base class of OmegaPhi3PiDecayer. */
  typedef FlatDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  OmegaPhi3PiDecayer class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<OmegaPhi3PiDecayer>
  : public ClassTraitsBase<OmegaPhi3PiDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::OmegaPhi3PiDecayer"; }
  /** Return the name of the shared library be loaded to get access to
   *  the OmegaPhi3PiDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "OmegaPhi3PiDecayer.so"; }
};

/** @endcond */

}

#endif /* THEPEG_OmegaPhi3PiDecayer_H */
