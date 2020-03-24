// -*- C++ -*-
//
// PtrTraits.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_PtrTraits_H
#define ThePEG_PtrTraits_H
// This is the declaration of the PtrTraits class.

namespace ThePEG {
namespace Pointer {

/**
 * PtrTraitsType is an empty non-polymorphic base class for all
 * PtrTraits classes.
 */
struct PtrTraitsType {};

/**
 * The PtrTraits class is used everywhere in ThePEG to
 * interface to the pointers which are handled. In particular, ThePEG
 * never uses new or delete but always
 * PtrTraits<P>::create and
 * PtrTraits<P>::destroy (to be precise the destroy method
 * is never used since all pointers are assumed to be reference
 * counted or in another way garbage collected). Also ThePEG
 * always uses dynamic_ptr_cast (rather than the standard
 * dynamic_cast) which in turn calls the
 * PtrTraits<P>::DynamicCast.
 *
 * In this file is also defined the specialized std::iterator_traits
 * for the reference counted pointers.
 *
 */
template <class T>
struct PtrTraits: public PtrTraitsType {};

/**
 * Specialization of the PtrTraits class for standard bare pointers.
 */
template <class T>
struct PtrTraits<T *>: public PtrTraitsType {

  /** Template argument typedef. */
  typedef T value_type;
  /** Template argument typedef. */
  typedef T & reference;
  /** Template argument typedef. */
  typedef const T & const_reference;
  /** Template argument typedef. */
  typedef T * pointer;
  /** Template argument typedef. */
  typedef T * const_pointer;

  /**
   * Return the bare pointer of the given pointer object.
   */
  static T * barePointer(T * p) { return p; }

  /**
   * Create an object and return a pointer to it.
   */
  static pointer create() { return new T; }

  /**
   * Create an copy of an object and return a pointer to it.
   */
  static pointer create(const_reference t) { return new T(t); }

  /**
   * Destroy the object pointed to.
   */
  static void destroy(pointer tp) { delete tp; }

  /**
   * Cast dynamically.
   */
  template <class R>
  static pointer DynamicCast(R * r) { return dynamic_cast<pointer>(r); }

  /**
   * Cast away constness.
   */
  static pointer ConstCast(const T * t) { return const_cast<pointer>(t); }

  /**
   * Cast from a basic pointer.
   */
  static pointer PtrCast(T * t) { return t; }

  /**
   * The bare pointer is not reference counted.
   */
  static const bool reference_counted = false;

};

/**
 * Specialization of the PtrTraits class for standard bare
 * const pointers.
 */
template <class T>
struct PtrTraits<const T *>: public PtrTraitsType {

  /** Template argument typedef. */
  typedef T value_type;
  /** Template argument typedef. */
  typedef T & reference;
  /** Template argument typedef. */
  typedef const T & const_reference;
  /** Template argument typedef. */
  typedef T * pointer;
  /** Template argument typedef. */
  typedef T * const_pointer;

  /**
   * Return the bare pointer of the given pointer object.
   */
  static const T * barePointer(const T * p) { return p; }

  /**
   * Create an object and return a pointer to it.
   */
  static pointer create() { return new T; }

  /**
   * Create an copy of an object and return a pointer to it.
   */
  static pointer create(const_reference t) { return new T(t); }

  /**
   * Destroy the object pointed to.
   */
  static void destroy(pointer tp) { delete tp; }

  /**
   * Cast dynamically.
   */
  template <class R>
  static const_pointer DynamicCast(const R * r) {
    return dynamic_cast<const_pointer>(r);
  }

  /**
   * Do not cast away constness.
   */
  static const_pointer ConstCast(const T * r) { return r; }

  /**
   * Cast from a basic pointer.
   */
  static const_pointer PtrCast(const T * t) { return t; }

  /**
   * The bare pointer is not reference counted.
   */
  static const bool reference_counted = false;

};

/**
 * Replacement for the standard dynamic_cast
 */
template <class T1, class T2>
T1 dynamic_ptr_cast(const T2 & t2) { return PtrTraits<T1>::DynamicCast(t2); }


/**
 * Replacement for the standard const_cast
 */
template <class T1, class T2>
T1 const_ptr_cast(const T2 & t2) { return PtrTraits<T1>::ConstCast(t2); }

/**
 * Simple interface to the PtrTraits<Ptr>::create()
 */
template <typename Ptr>
inline Ptr ptr_new() { return PtrTraits<Ptr>::create(); }

/**
 * Simple interface to the PtrTraits<Ptr>::create()
 */
template <typename Ptr>
inline Ptr ptr_new(typename PtrTraits<Ptr>::const_reference t) {
  return PtrTraits<Ptr>::create(t);
}

/**
 * Simple interface to the PtrTraits<Ptr>::create()
 */
template <typename T>
inline typename Ptr<T>::pointer new_ptr() {
  return PtrTraits< typename Ptr<T>::pointer >::create();
}

/**
 * Simple interface to the PtrTraits<Ptr>::create()
 */
template <typename T>
inline typename Ptr<T>::pointer new_ptr(const T & t) {
  return PtrTraits< typename Ptr<T>::pointer >::create(t);
}

/**
 * Simple interface to the PtrTraits<Ptr>::PtrCast()
 */
template <typename TPtr, typename T>
inline TPtr ptr_cast(T * t) {
  return PtrTraits<TPtr>::PtrCast(t);
}

/**
 * Simple interface to the PtrTraits<Ptr>::PtrCast()
 */
template <typename TPtr, typename T>
inline TPtr ptr_cast_const(const T * t) {
  return PtrTraits<TPtr>::PtrCast(const_cast<T*>(t));
}


}
}

#endif /* ThePEG_PtrTraitsH */
