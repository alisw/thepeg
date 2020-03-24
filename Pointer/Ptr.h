// -*- C++ -*-
//
// Ptr.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Ptr_H
#define ThePEG_Ptr_H
// This is the declaration of the Ptr class.

#include "RCPtr.fh"

namespace ThePEG {

/** The namespace for the reference counted pointer classes */
namespace Pointer {

/**
 * Ptr is a templated class to provide typedefs for pointers types
 * ThePEG should use for a given type. If you want to use ThePEG with
 * another kind of (smart) pointers than those provided you have to
 * provide an alternative ThePEG.h file which includes an alternative
 * Ptr class. It may also be possible to specialize this Ptr class for
 * using different pointer classes for different classes.
 *
 * @see RCPtr
 * @see ConstRCPtr
 * @see TransientRCPtr
 * @see TransientConstRCPtr
 * 
 */
template <typename T>
struct Ptr {

  /** Template argument typedef. */
  typedef RCPtr<T> pointer;
  /** Template argument typedef. */
  typedef ConstRCPtr<T> const_pointer;
  /** Template argument typedef. */
  typedef TransientRCPtr<T> transient_pointer;
  /** Template argument typedef. */
  typedef TransientConstRCPtr<T> transient_const_pointer;
  /** Template argument typedef. */
  typedef pointer ptr;
  /** Template argument typedef. */
  typedef const_pointer cptr;
  /** Template argument typedef. */
  typedef transient_pointer tptr;
  /** Template argument typedef. */
  typedef transient_const_pointer tcptr;
  /** Template argument typedef. */
  typedef pointer p;
  /** Template argument typedef. */
  typedef const_pointer cp;
  /** Template argument typedef. */
  typedef transient_pointer tp;
  /** Template argument typedef. */
  typedef transient_const_pointer tcp;

};

}
}

#endif /* ThePEG_Ptr_H */
