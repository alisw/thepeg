// -*- C++ -*-
//
// StandardRandom.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_StandardRandom_H
#define ThePEG_StandardRandom_H
// This is the declaration of the StandardRandom class.

#include "RandomGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace ThePEG {

/**
 * StandardRandom inherits from the RandomGenerator class and is an
 * interface to the CLHEP::JamesRandom engine class.
 *
 * @see \ref StandardRandomInterfaces "The interfaces"
 * defined for StandardRandom.
 */
class StandardRandom: public RandomGenerator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  StandardRandom() : u() { if ( theSeed != 0 ) setSeed(theSeed); }
  //@}

public:

  /**
   * Reset the underlying random algorithm with the given seed. If the
   * \a seed is set to -1 a standard seed will be used.
   */
  virtual void setSeed(long seed);

protected:

  /**
   * Fill the cache with random numbers.
   */
  virtual void fill();

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
   * The internal state vector.
   */
  array<double,97> u;

  /**
   * Parameter for the internal state.
   */
  double c;

  /**
   * Parameter for the internal state.
   */
  double cd;

  /**
   * Parameter for the internal state.
   */
  double cm;

  /**
   * Index for the internal state.
   */
  int i97;

  /**
   * Index for the internal state.
   */
  int j97;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<StandardRandom> initStandardRandom;

  /**
   *  Private and non-existent assignment operator.
   */
  StandardRandom & operator=(const StandardRandom &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of StandardRandom. */
template <>
struct BaseClassTrait<StandardRandom,1>: public ClassTraitsType {
  /** Typedef of the first base class of StandardRandom. */
  typedef RandomGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  StandardRandom class. */
template <>
struct ClassTraits<StandardRandom>: public ClassTraitsBase<StandardRandom> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::StandardRandom"; }
};

/** @endcond */

}

#endif /* ThePEG_StandardRandom_H */
