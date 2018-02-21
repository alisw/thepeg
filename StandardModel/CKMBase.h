// -*- C++ -*-
//
// CKMBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_CKMBase_H
#define ThePEG_CKMBase_H
// This is the declaration of the CKMBase class.

#include "ThePEG/Interface/Interfaced.h"
#include "StandardModelBase.fh"

namespace ThePEG {

/**
 * CKMBase is an abstract base classused by the StandardModelBase to
 * implement the Cabibbo-Kobayashi-Mascawa mixing matrix. Concrete
 * sub-classes must implement the getMatrix(unsigned int) function.
 *
 * @see \ref CKMBaseInterfaces "The interfaces"
 * defined for CKMBase.
 * @see StandardModelBase
 */
class CKMBase: public Interfaced {

public:

  /**
   * Return the matrix of squared CKM matrix elements. The returned
   * matrix should be for \a nf families. This function must be
   * overridden by sub-classes.
   */
  virtual vector< vector<double> >  getMatrix(unsigned int nf) const = 0;

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

private:

  /**
   * Describe an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<CKMBase> initCKMBase;

  /**
   *  Private and non-existent assignment operator.
   */
  CKMBase & operator=(const CKMBase &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of CKMBase. */
template <>
struct BaseClassTrait<CKMBase,1>: public ClassTraitsType {
  /** Typedef of the first base class of CKMBase. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  CKMBase class. */
template <>
struct ClassTraits<CKMBase>: public ClassTraitsBase<CKMBase> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::CKMBase"; }
};

/** @endcond */

}

#endif /* ThePEG_CKMBase_H */
