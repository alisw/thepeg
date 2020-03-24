// -*- C++ -*-
//
// ReweightConstant.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_ReweightConstant_H
#define THEPEG_ReweightConstant_H
//
// This is the declaration of the ReweightConstant class.
//

#include "ThePEG/MatrixElement/ReweightBase.h"

namespace ThePEG {

/**
 * The ReweightConstant class is a simple ReweightBase sub-class which
 * simply reweight an event with a constant.
 *
 * @see \ref ReweightConstantInterfaces "The interfaces"
 * defined for ReweightConstant.
 */
class ReweightConstant: public ReweightBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ReweightConstant() : C(1.0) {}
  //@}

public:

  /**
   * Return the wieght for the kinematical configuation provided by
   * the assigned XComb object (in the LastXCombInfo base class).
   */
  virtual double weight() const;

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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
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
   * The constant to reweight with.
   */
  double C;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ReweightConstant> initReweightConstant;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ReweightConstant & operator=(const ReweightConstant &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ReweightConstant. */
template <>
struct BaseClassTrait<ReweightConstant,1> {
  /** Typedef of the first base class of ReweightConstant. */
  typedef ReweightBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ReweightConstant class and the shared object where it is defined. */
template <>
struct ClassTraits<ReweightConstant>
  : public ClassTraitsBase<ReweightConstant> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::ReweightConstant"; }
  /** Return the name of the shared library be loaded to get
   *  access to the ReweightConstant class and every other class it uses
   *  (except the base class). */
  static string library() { return "ReweightConstant.so"; }
};

/** @endcond */

}

#endif /* THEPEG_ReweightConstant_H */
