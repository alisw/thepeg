// -*- C++ -*-
//
// WeakToHadronsDecayer.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_WeakToHadronsDecayer_H
#define THEPEG_WeakToHadronsDecayer_H
// This is the declaration of the WeakToHadronsDecayer class.

#include "ThePEG/PDT/QuarksToHadronsDecayer.h"

namespace ThePEG {

/**
 * The WeakToHadronsDecayer class inherits from QuarksToHadronsDecayer
 * and can performs weak decays of taus and charmed and bottom
 * hadrons. The intermediate W can either decay leptonically in which
 * case standard V-A matrix element is used, or it can decay into
 * quarks in which case the conversion into quarks is performed as for
 * the QuarkToHadronsDecayer base class. In both cases the W decay
 * products should be specified first. The spectator system can either
 * be specified in terms of hadrons or in terms of quarks which will
 * be collapsed into a single hadron.
 *
 * @see \ref WeakToHadronsDecayerInterfaces "The interfaces"
 * defined for WeakToHadronsDecayer.
 */
class WeakToHadronsDecayer: public QuarksToHadronsDecayer {

public:

  /** @name Virtual functions required by the Decayer and
      QuarksToHadronsDecayer classes.
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
   * Called by QuarksToHadronsDecayer::distribute() to reweight the
   * default flat phase spece. Can be overridden by sub-classes and
   * should return a number between 0 and 1. This version returns 1.
   */
  virtual double reweight(const Particle & parent,
			  const PVector & children) const;

  /**
   * Produce \a Nh hadrons from the specified \a quarks. The last
   * quark is considered to be a spectator quark.
   */
  virtual PVector getHadrons(int Nh, tcPDVector quarks) const;
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
  static ClassDescription<WeakToHadronsDecayer> initWeakToHadronsDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  WeakToHadronsDecayer & operator=(const WeakToHadronsDecayer &) = delete;

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of WeakToHadronsDecayer. */
template <>
struct BaseClassTrait<WeakToHadronsDecayer,1>: public ClassTraitsType {
  /** Typedef of the first base class of WeakToHadronsDecayer. */
  typedef QuarksToHadronsDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  WeakToHadronsDecayer class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<WeakToHadronsDecayer>
  : public ClassTraitsBase<WeakToHadronsDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::WeakToHadronsDecayer"; }
  /** Return the name of the shared library be loaded to get access to
   *  the WeakToHadronsDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "WeakToHadronsDecayer.so"; }
};

/** @endcond */

}

#endif /* THEPEG_WeakToHadronsDecayer_H */
