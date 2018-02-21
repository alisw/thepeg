// -*- C++ -*-
//
// Complex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Complex_H
#define ThePEG_Complex_H
//
// This is file wraps the standard complex header and makes some
// convenient typedefs in the ThePEG namespace.
//

#include <complex>

namespace ThePEG {

using std::complex;

/** ThePEG code should use Complex for all complex scalars */
typedef std::complex<double> Complex;

/** @cond TRAITSPECIALIZATIONS */

template <typename T, typename U>
struct BinaryOpTraits;


template <typename T, typename U>
struct BinaryOpTraits<complex<T>, U> {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef complex<typename BinaryOpTraits<T,U>::MulT> MulT;
  /** The type resulting from division of one template type with
      another. */
  typedef complex<typename BinaryOpTraits<T,U>::DivT> DivT;
  
};

template <typename T, typename U>
struct BinaryOpTraits<T, complex<U> > {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef complex<typename BinaryOpTraits<T,U>::MulT> MulT;
  /** The type resulting from division of one template type with
      another. */
  typedef complex<typename BinaryOpTraits<T,U>::DivT> DivT;
  
};

template <typename T, typename U>
struct BinaryOpTraits<complex<T>, complex<U> > {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef complex<typename BinaryOpTraits<T,U>::MulT> MulT;
  /** The type resulting from division of one template type with
      another. */
  typedef complex<typename BinaryOpTraits<T,U>::DivT> DivT;
  
};

template <typename T>
struct BinaryOpTraits<complex<T>, complex<T> > {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef complex<typename BinaryOpTraits<T,T>::MulT> MulT;
  /** The type resulting from division of one template type with
      another. */
  typedef complex<typename BinaryOpTraits<T,T>::DivT> DivT;

  /** @endcond */

};



}

#endif /* ThePEG_Complex_H */

