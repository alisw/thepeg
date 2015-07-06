// -*- C++ -*-
//
// ThePEGStrategy.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ThePEGStrategy_H
#define ThePEG_ThePEGStrategy_H
// This is the declaration of the ThePEGStrategy class.

#include "ThePEG/Repository/Strategy.h"

namespace ThePEG {

/**
 * The ThePEGStrategy class is a sub-class of the Strategy class,
 * simply implementing the correct citation for ThePEG in the
 * ClassDocumentation interface.
 *
 * @see Strategy
 * 
 */
class ThePEGStrategy: public Strategy {

public:

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
   * Describe concrete class without persistent data.
   */
  static NoPIOClassDescription<ThePEGStrategy> initThePEGStrategy;

  /**
   *  Private and non-existent assignment operator.
   */
  ThePEGStrategy & operator=(const ThePEGStrategy &);

};


/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of ThePEGStrategy. */
template <>
struct BaseClassTrait<ThePEGStrategy,1>: public ClassTraitsType {
  /** Typedef of the first base class of ThePEGStrategy. */
  typedef Strategy NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  ThePEGStrategy class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<ThePEGStrategy>: public ClassTraitsBase<ThePEGStrategy> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::ThePEGStrategy"; }
  /** Return the name of the shared library be loaded to get access to
   *  the ThePEGStrategy class and every other class it uses
   *  (except the base class). */
  static string library() { return "ThePEGStrategy.so"; }
};

/** @endcond */

}

#endif /* ThePEG_ThePEGStrategy_H */
