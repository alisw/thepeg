// -*- C++ -*-
//
// SimpleZGenerator.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_SimpleZGenerator_H
#define THEPEG_SimpleZGenerator_H
// This is the declaration of the SimpleZGenerator class.

#include "ThePEG/Handlers/ZGenerator.h"

namespace ThePEG {

/**
 * SimpleZGenerator is a very simple concrete subclass of
 * ZGenerator. It implements a naive unphysical model to generate the
 * momentum fraction, \f$z\f$, taken by hadrons produced in a hadronization
 * scenario. It should only be used for testing purposes.
 *
 * @see \ref SimpleZGeneratorInterfaces "The interfaces"
 * defined for SimpleZGenerator.
 */
class SimpleZGenerator: public ZGenerator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Destructor.
   */
  virtual ~SimpleZGenerator();
  //@}

public:

  /** @name Virtual functions mandated by the ZGenerator base class. */
  //@{
  /**
   * Return the momentum fraction. Assume that an initial
   * (anti-)(di-)quark \a q1 produces a hadron and leaves behind
   * another (anti-)(di-)quark \a q2. The hadron is assumed to have a
   * squared transverse mass, \a mT2, w.r.t. the initial quark
   * direction.
   * @return the energy fraction distributed as \f$\sqrt{z}\f$ (or
   * \f$1-\sqrt{z}\f$) if \a q1 (or \a q2) is a diquark. Otherwise a
   * flat distribution is used.
   */
  virtual double generate(cPDPtr q1, cPDPtr q2, Energy2 mT2 ) const;
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
  static ClassDescription<SimpleZGenerator> initSimpleZGenerator;

  /**
   * Private and non-existent assignment operator.
   */
  SimpleZGenerator & operator=(const SimpleZGenerator &) = delete;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * SimpleZGenerator.
 */
template <>
struct BaseClassTrait<SimpleZGenerator,1>: public ClassTraitsType {
  /** Typedef of the base class of SimpleZGenerator. */
  typedef ZGenerator NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * SimpleZGenerator class and the shared object where it is defined.
 */
template <>
struct ClassTraits<SimpleZGenerator>
  : public ClassTraitsBase<SimpleZGenerator> {
  /** Return the class name.  */
  static string className() { return "ThePEG::SimpleZGenerator"; }
  /**
   * Return the name of the shared library to be loaded to get access
   * to the SimpleZGenerator class and every other class it uses
   * (except the base class).
   */
  static string library() { return "SimpleZGenerator.so"; }

};

/** @endcond */

}

#endif /* THEPEG_SimpleZGenerator_H */
