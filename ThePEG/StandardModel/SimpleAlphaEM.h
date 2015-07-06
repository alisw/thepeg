// -*- C++ -*-
//
// SimpleAlphaEM.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_SimpleAlphaEM_H
#define ThePEG_SimpleAlphaEM_H
// This is the declaration of the SimpleAlphaEM class.

#include "AlphaEMBase.h"

namespace ThePEG {

/**
 * SimpleAlphaEM inherits from AlphaEMBase and implements a simple
 * running of the electromagnetic coupling as parameterized by
 * H.~Buckhardt et al.
 *
 * @see \ref SimpleAlphaEMInterfaces "The interfaces"
 * defined for SimpleAlphaEM.
 */
class SimpleAlphaEM: public AlphaEMBase {

public:

  /**
   * The \f$\alpha_{EM}\f$. Return the value of the coupling at a
   * given \a scale using the given standard model object, \a sm.
   */
  virtual double value(Energy2 scale, const StandardModelBase &) const;

  /**
   * Return the number of loops contributing to
   * the running this coupling.
   */
  virtual unsigned int nloops () const { return 1; }

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
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<SimpleAlphaEM> initSimpleAlphaEM;

  /**
   *  Private and non-existent assignment operator.
   */
  SimpleAlphaEM & operator=(const SimpleAlphaEM &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of SimpleAlphaEM. */
template <>
struct BaseClassTrait<SimpleAlphaEM,1>: public ClassTraitsType {
  /** Typedef of the first base class of SimpleAlphaEM. */
  typedef AlphaEMBase NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  SimpleAlphaEM class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<SimpleAlphaEM>: public ClassTraitsBase<SimpleAlphaEM> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::SimpleAlphaEM"; }
  /** Return the name of the shared library be loaded to get access to
   *  the SimpleAlphaEM class and every other class it uses
   *  (except the base class). */
  static string library() { return "SimpleAlphaEM.so"; }
};

/** @endcond */

}

#endif /* ThePEG_SimpleAlphaEM_H */
