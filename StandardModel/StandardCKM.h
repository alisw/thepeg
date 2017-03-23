// -*- C++ -*-
//
// StandardCKM.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_StandardCKM_H
#define ThePEG_StandardCKM_H
// This is the declaration of the StandardCKM class.

#include "CKMBase.h"

namespace ThePEG {

/**
 * StandardCKM inherits from CKMBase and implements the standard
 * parameterization of the CKM matrix in terms of three angles and a
 * phase.
 *
 * @see \ref StandardCKMInterfaces "The interfaces"
 * defined for StandardCKM.
 */
class StandardCKM: public CKMBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  StandardCKM() : theta12(0.222357), theta13(0.0003150), 
		  theta23(0.039009), delta(1.35819) {}
  //@}

public:

  /**
   * Return the matrix of squared CKM matrix elements. The returned
   * matrix should be for \a nf families.
   */
  virtual vector< vector<double> >  getMatrix(unsigned int nFamilies) const;

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
   * The \f$\theta_{12}\f$ angle.
   */
  double theta12;

  /**
   * The \f$\theta_{13}\f$ angle.
   */
  double theta13;

  /**
   * The \f$\theta_{23}\f$ angle.
   */
  double theta23;

  /**
   * The \f$\delta\f$ angle describing the phase.
   */
  double delta;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<StandardCKM> initStandardCKM;

  /**
   *  Private and non-existent assignment operator.
   */
  StandardCKM & operator=(const StandardCKM &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of StandardCKM. */
template <>
struct BaseClassTrait<StandardCKM,1>: public ClassTraitsType {
  /** Typedef of the first base class of StandardCKM. */
  typedef CKMBase NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  StandardCKM class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<StandardCKM>: public ClassTraitsBase<StandardCKM> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::StandardCKM"; }
  /** Return the name of the shared library be loaded to get access to
   *  the StandardCKM class and every other class it uses
   *  (except the base class). */
  static string library() { return "StandardCKM.so"; }
};

/** @endcond */

}

#endif /* ThePEG_StandardCKM_H */
