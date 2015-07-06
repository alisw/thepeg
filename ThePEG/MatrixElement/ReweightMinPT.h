// -*- C++ -*-
//
// ReweightMinPT.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ReweightMinPT_H
#define ThePEG_ReweightMinPT_H
// This is the declaration of the ReweightMinPT class.

#include "ThePEG/MatrixElement/ReweightBase.h"

namespace ThePEG {

/**
 * The ReweightMinPT class reweights matrix elements with the minimum
 * of the transverse momenta of the outgoing partons to some power.
 *
 * @see ReweightBase
 * 
 */
class ReweightMinPT: public ReweightBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  ReweightMinPT()
    : power(4.0), scale(50.0*GeV) {}
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
   * The weight is the minimum pt/scale to a \a power.
   */
  double power;

  /**
   * The weight is the minimum pt/\a scale to a power.
   */
  Energy scale;

private:

  /**
   * Describe a concrete base class with persistent data.
   */
  static ClassDescription<ReweightMinPT> initReweightMinPT;

  /**
   *  Private and non-existent assignment operator.
   */
  ReweightMinPT & operator=(const ReweightMinPT &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ReweightMinPT. */
template <>
struct BaseClassTrait<ReweightMinPT,1>: public ClassTraitsType {
  /** Typedef of the first base class of ReweightMinPT. */
  typedef ReweightBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ReweightMinPT class and the shared object where it is defined. */
template <>
struct ClassTraits<ReweightMinPT>: public ClassTraitsBase<ReweightMinPT> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::ReweightMinPT"; }
  /** Return the name of the shared library be loaded to get
   *  access to the ReweightMinPT class and every other class it uses
   *  (except the base class). */
  static string library() { return "ReweightMinPT.so"; }
};

/** @endcond */

}

#endif /* ThePEG_ReweightMinPT_H */
