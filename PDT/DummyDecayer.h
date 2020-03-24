// -*- C++ -*-
//
// DummyDecayer.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_DummyDecayer_H
#define THEPEG_DummyDecayer_H
// This is the declaration of the DummyDecayer class.

#include "ThePEG/PDT/Decayer.h"

namespace ThePEG {

/**
 * DummyDecayer inherits from Decayer and is a dummy decayer class to
 * be used for symbolic decay channels. If it for some reason is
 * called to perform a decay, it will throw a std::logic_error.
 *
 * @see \ref DummyDecayerInterfaces "The interfaces"
 * defined for DummyDecayer.
 */
class DummyDecayer: public Decayer {

public:

  /** @name Virtual functions required by the Decayer class.
   */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true always.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * Perform a decay for a given DecayMode and a given Particle
   * instance. Will throw std::logic_error if called.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;
 
  /**
   * Will always return zero, since no decay can ever be performed
   * with this decayer.
   */
  virtual double brat(const DecayMode &,
		      const ParticleData &, double) const;
  /**
   * Will always return zero, since no decay can ever be performed
   * with this decayer.
   */
  virtual double brat(const DecayMode &, const Particle &, double) const;
  //@}

public:

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
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<DummyDecayer> initDummyDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  DummyDecayer & operator=(const DummyDecayer &) = delete;

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of DummyDecayer. */
template <>
struct BaseClassTrait<DummyDecayer,1>: public ClassTraitsType {
  /** Typedef of the first base class of DummyDecayer. */
  typedef Decayer NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  DummyDecayer class. */
template <>
struct ClassTraits<DummyDecayer>
  : public ClassTraitsBase<DummyDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::DummyDecayer"; }

};

/** @endcond */

}

#endif /* THEPEG_DummyDecayer_H */
