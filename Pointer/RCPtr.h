// -*- C++ -*-
//
// RCPtr.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_RCPtr_H
#define ThePEG_RCPtr_H
// This is the declaration of the RCPtrBase,


#include "ReferenceCounted.h"
#include "RCPtr.fh"
#include "PtrTraits.h"

namespace ThePEG {
namespace Pointer {

/**
 * RCPtrBase is the base class of RCPtr and ConstRCPtr which are
 * reference counted (smart) pointers. The RCPtrBase class communicates with the
 * ReferenceCounted object which must be the base class of
 * all objects pointed to and which keeps count of the pointers
 * pointing to an object.
 *
 * @see ReferenceCounted
 */
class RCPtrBase {

  /** Get counter type from ReferenceCounted class. */
  typedef ReferenceCounted::CounterType CounterType;

protected:

  /**
   * Increment the counter of a reference counted object.
   */
  void increment(const ReferenceCounted * rcp) {
    if ( rcp ) rcp->incrementReferenceCount();
  }
  /**
   * Decrement the counter of a reference counted object.
   */
  bool release(const ReferenceCounted * rcp) {
    return rcp && rcp->decrementReferenceCount();
  }

};

/**
 * RCPtr is a reference counted (smart) pointer. Objects created using
 * the RCPtr::create() methods will continue living until no RCPtr or
 * ConstRCPtr are pointing to it, at which point it will be deleted.
 *
 * @see ReferenceCounted
 */
template <typename T>
class RCPtr: public RCPtrBase {

public:

  /** Template argument typedef. */
  typedef void iterator_category;
  /** Template argument typedef. */
  typedef void difference_type;
  /** Template argument typedef. */
  typedef T * pointer;
  /** Template argument typedef. */
  typedef const T * const_pointer;
  /** Template argument typedef. */
  typedef T & reference;
  /** Template argument typedef. */
  typedef const T & const_reference;
  /** Template argument typedef. */
  typedef T value_type;

public:

  /** <code></code> */
  /**
   * Constructor for null pointer.
   */
  RCPtr() : ptr(nullptr) {}

  /**
   * Constructor for nullptr.
   */
  RCPtr( decltype(nullptr) ) : ptr(nullptr) {}

  /**
   * Copy constructor.
   */
  RCPtr(const RCPtr & p) : ptr(p.ptr) { increment(); }

  /**
   * Copy constructor for class UPtr which has operator-> defined
   * resulting in a value implicitly convertible to T *.
   */
  template <typename UPtr>
  RCPtr(const UPtr & u) 
    : ptr(PtrTraits<UPtr>::barePointer(u)) { increment(); }

  /**
   * Construction from real pointer.
   */
  explicit RCPtr(pointer p) : ptr(p) { increment(); }

  /**
   * Destructor. Deletes the object pointed to if this is the last
   * pointer to it.
   */
  ~RCPtr() { release(); }

  /**
   * Allocate and construct an object of class T and return a RCPtr to
   * it.
   */
  static RCPtr Create() {
    RCPtr<T> p;
    return p.create();
  }

  /**
   * Allocate and copy-construct an object of class T and return a
   * RCPtr to it.
   */
  static RCPtr Create(const_reference t) {
    RCPtr<T> p;
    return p.create(t);
  }


  /**
   * Allocate and construct an object of class T and point to it,
   * possibly deleting the object pointed to before.
   */
  RCPtr & create() {
    release();
    ptr = new T;
    //  increment(); // ReferenceCounted() constructor starts at 1
    return *this;
  }

  /**
   * Allocate and copy-construct an object of class T and point to it,
   * possibly deleting the object pointed to before.
   */
  RCPtr & create(const_reference t) {
    release();
    ptr = new T(t);
    //  increment(); // ReferenceCounted() constructor starts at 1
    return *this;
  }

  /**
   * Assignment.
   */
  RCPtr & operator=(const RCPtr & p) {
    if ( ptr == p.ptr ) return *this;
    release();
    ptr = p.ptr;
    increment();
    return *this;
  }

  /**
   * Assignment from class UPtr which has operator-> defined resulting
   * in a value implicitly convertible to T *.
   */
  template <typename UPtr>
  RCPtr & operator=(const UPtr & u) {
    if ( ptr == PtrTraits<UPtr>::barePointer(u) ) return *this;
    release();
    ptr = PtrTraits<UPtr>::barePointer(u);
    increment();
    return *this;
  }

  /**
   * Assignment from class UPtr which has operator-> defined resulting
   * in a value which can be cast dynamically to T *.
   */
  template <typename UPtr>
  RCPtr & assignDynamic(const UPtr & u) {
    pointer up = dynamic_cast<pointer>(PtrTraits<UPtr>::barePointer(u));
    if ( ptr == up ) return *this;
    release();
    ptr = up;
    increment();
    return *this;
  }

  /**
   * Assignment from class UPtr which has operator-> defined resulting
   * in a value which can be const_cast'ed to T *.
   */
  template <typename UPtr>
  RCPtr & assignConst(const UPtr & u) {
    pointer up = const_cast<pointer>(PtrTraits<UPtr>::barePointer(u));
    if ( ptr == up ) return *this;
    release();
    ptr = up;
    increment();
    return *this;
  }

  /**
   * Make p point to the object pointed to by this and vice versa.  
   */
  void swap(RCPtr & p) {
    const pointer tmp = ptr;
    ptr = p.ptr;
    p.ptr = tmp;
    //  std::swap(ptr, p.ptr);
  }

  /**
   * Test for equality of the underlying pointers.
   */
  bool operator==(const RCPtr & p) const { return ptr == p.ptr; }

  /**
   * Test for inequality of the underlying pointers.
   */
  bool operator!=(const RCPtr & p) const { return ptr != p.ptr; }

  /**
   * Test for equality of the underlying pointers.
   */
  bool operator==(const_pointer p) const { return ptr == p; }

  /**
   * Test for inequality of the underlying pointers.
   */
  bool operator!=(const_pointer p) const { return ptr != p; }

  /**
   * Test for equality of the underlying pointers.
   */
  template <typename UPtr>
  bool operator==(const UPtr & u) const {
    return ptr == PtrTraits<UPtr>::barePointer(u);
  }

  /**
   * Test for inequality of the underlying pointers.
   */
  template <typename UPtr>
  bool operator!=(const UPtr & u) const {
    return ptr != PtrTraits<UPtr>::barePointer(u);
  }

  /**
   * Test for ordering.
   */
  bool operator<(const RCPtr & p) const {
    return ( ptr && p.ptr && ptr->uniqueId != p.ptr->uniqueId ) ?
      ptr->uniqueId < p.ptr->uniqueId : ptr < p.ptr;
  }

  /**
   * Test for ordering.
   */
  bool operator<(const_pointer p) const {
    return ( ptr && p && ptr->uniqueId != p->uniqueId ) ?
      ptr->uniqueId < p->uniqueId : ptr < p;
  }

  /**
   * Returns true if the underlying pointer is null.
   */
  bool operator!() const { return !ptr; }

  /**
   * Returns the underlying pointer.
   */
  operator T * () const { return ptr; }

  /**
   * Member access.
   */
  pointer operator->() const { return ptr; }

  /**
   * Dereferencing.
   */
  reference operator*() const { return *ptr; }

private:

  /**
   * Increment the counter of the object pointed to.
   */
  void increment() { RCPtrBase::increment(ptr); }

  /**
   * Stop pointing to the current object and delete it if this was the
   * last pointer to it.
   */
  void release() { if ( RCPtrBase::release(ptr) )  delete ptr; }
   
  /**
   * The actual pointer.
   */
  pointer ptr;

};

/**
 * ConstRCPtr is a reference counted (smart) const pointer. Objects
 * created using the RCPtr::create() methods will continue living
 * until no RCPtr or ConstRCPtr are pointing to it, at which point it
 * will be deleted.
 *
 * @see ReferenceCounted
 */
template <typename T>
class ConstRCPtr : public RCPtrBase {

public:

  /** Template argument typedef. */
  typedef void iterator_category;
  /** Template argument typedef. */
  typedef void difference_type;
  /** Template argument typedef. */
  typedef T * pointer;
  /** Template argument typedef. */
  typedef const T * const_pointer;
  /** Template argument typedef. */
  typedef T & reference;
  /** Template argument typedef. */
  typedef const T & const_reference;
  /** Template argument typedef. */
  typedef T value_type;

public:

  /**
   * Constructor for null pointer.
   */
  ConstRCPtr() : ptr(nullptr) {}

  /**
   * Constructor for nullptr.
   */
  ConstRCPtr( decltype(nullptr) ) : ptr(nullptr) {}

  /**
   * Copy constructor.
   */
  ConstRCPtr(const ConstRCPtr & p) : ptr(p.ptr) { increment(); }

  /**
   * Copyconstructor for class UPtr which has operator-> defined
   * resulting in a value implicitly convertible to const T *.
   */
  template <typename UPtr>
  ConstRCPtr(const UPtr & u) : ptr(PtrTraits<UPtr>::barePointer(u)) { increment(); }

  /**
   * Construction from real pointer.
   */
  explicit ConstRCPtr(const_pointer p) : ptr(p) { increment(); }

  /**
   * Destructor. Deletes the object pointed to if this is the last
   * pointer to it.
   */
  ~ConstRCPtr() { release(); }

  /**
   * Assignment.
   */
  ConstRCPtr & operator=(const ConstRCPtr & p) {
    if ( ptr == p.ptr ) return *this;
    release();
    ptr = p.ptr;
    increment();
    return *this;
  }
  
  /**
   * Assignment from class UPtr which has operator-> defined resulting
   * in a value implicitly convertible to const T *.
   */
  template <typename UPtr>
  ConstRCPtr & operator=(const UPtr & u) {
    if ( ptr == PtrTraits<UPtr>::barePointer(u) ) return *this;
    release();
    ptr = PtrTraits<UPtr>::barePointer(u);
    increment();
    return *this;
  }

  /**
   * Assignment from class UPtr which has operator-> defined resulting
   * in a value which can be cast dynamically to const T *.
   */
  template <typename UPtr>
  ConstRCPtr & assignDynamic(const UPtr & u) {
    const_pointer up =
      dynamic_cast<const_pointer>(PtrTraits<UPtr>::barePointer(u));
    if ( ptr == up ) return *this;
    release();
    ptr = up;
    increment();
    return *this;
  }

  /**
   * Make p point to the object pointed to by this and vice versa.  
   */
  void swap(ConstRCPtr & p) {
    const const_pointer tmp = ptr;
    ptr = p.ptr;
    p.ptr = tmp;
    //  std::swap(ptr, p.ptr);
  }  

  /**
   * Test for equality of the underlying pointers.
   */
  bool operator==(const ConstRCPtr & p) const { return ptr == p.ptr; }

  /**
   * Test for inequality of the underlying pointers.
   */
  bool operator!=(const ConstRCPtr & p) const { return ptr != p.ptr; }

  /**
   * Test for equality of the underlying pointers.
   */
  bool operator==(const_pointer p) const { return ptr == p; }

  /**
   * Test for inequality of the underlying pointers.
   */
  bool operator!=(const_pointer p) const { return ptr != p; }

  /**
   * Test for equality of the underlying pointers.
   */
  template <typename UPtr>
    bool operator==(const UPtr & u) const { return ptr == PtrTraits<UPtr>::barePointer(u); }

  /**
   * Test for inequality of the underlying pointers.
   */
  template <typename UPtr>
  bool operator!=(const UPtr & u) const { return ptr != PtrTraits<UPtr>::barePointer(u); }

  /**
   * Test for ordering.
   */
  bool operator<(const ConstRCPtr & p) const {
    return ( ptr && p.ptr && ptr->uniqueId != p.ptr->uniqueId ) ?
      ptr->uniqueId < p.ptr->uniqueId : ptr < p.ptr;
  }

  /**
   * Test for ordering.
   */
  bool operator<(const_pointer p) const {
    return ( ptr && p && ptr->uniqueId != p->uniqueId ) ?
      ptr->uniqueId < p->uniqueId : ptr < p;
  }

  /**
   * Returns true if the underlying pointer is null.
   */
  bool operator!() const { return !ptr; }

  /**
   * Returns the underlying pointer.
   */
  operator const T * () const { return ptr; }

  /**
   * Member access.
   */
  const_pointer operator->() const { return ptr; }

  /**
   * Dereferencing.
   */
  const_reference operator*() const { return *ptr; }

private:

  /**
   * Increment the counter of the object pointed to.
   */
  void increment() { RCPtrBase::increment(ptr); }

  /**
   * Stop pointing to the current object and delete it if this was the
   * last pointer to it.
   */
  void release() { if ( RCPtrBase::release(ptr) ) delete ptr; }
   
  /**
   * The actual pointer.
   */
  const_pointer ptr;

};

/**
 * TransientRCPtr is a simple wrapper around a bare pointer which can
 * be assigned to and from an RCPtr and ConstRCPtr without problem.
 *
 * @see RCPtr
 * @see ConstRCPtr
 */
template <typename T>
class TransientRCPtr {

public:

  /** Template argument typedef. */
  typedef void iterator_category;
  /** Template argument typedef. */
  typedef void difference_type;
  /** Template argument typedef. */
  typedef T * pointer;
  /** Template argument typedef. */
  typedef const T * const_pointer;
  /** Template argument typedef. */
  typedef T & reference;
  /** Template argument typedef. */
  typedef const T & const_reference;
  /** Template argument typedef. */
  typedef T value_type;

public:

  /**
   * Constructor for null pointer.
   */
  TransientRCPtr() : ptr(nullptr) {}

  /**
   * Constructor for nullptr.
   */
  TransientRCPtr( decltype(nullptr) ) : ptr(nullptr) {}

  /**
   * Copy constructor.
   */
  TransientRCPtr(const TransientRCPtr & p) : ptr(p.ptr) {}

  /**
   * Copyconstructor for class UPtr which has operator-> defined
   * resulting in a value implicitly convertible to T *.
   */
  template <typename UPtr>
  TransientRCPtr(const UPtr & u) : ptr(PtrTraits<UPtr>::barePointer(u)) {}

  /**
   * Construction from real pointer.
   */
  explicit TransientRCPtr(pointer p) : ptr(p) {}

  /**
   * Destructor. Does not delete the object pointed to.
   */
  ~TransientRCPtr() {}

  /**
   * Assignment.
   */
  TransientRCPtr & operator=(const TransientRCPtr & p) {
    ptr = p.ptr;
    return *this;
  }

  /**
   * Assignment from class UPtr which has operator-> defined resulting
   * in a value implicitly convertible to T *.
   */
  template <typename UPtr>
  TransientRCPtr & operator=(const UPtr & u) {
    ptr = PtrTraits<UPtr>::barePointer(u);
    return *this;
  }

  /**
   * Assignment from class UPtr which has operator-> defined resulting
   * in a value which can be cast dynamically to T *.
   */
  template <typename UPtr>
  TransientRCPtr & assignDynamic(const UPtr & u) {
    ptr = dynamic_cast<pointer>(PtrTraits<UPtr>::barePointer(u));
    return *this;
  }

  /**
   * Assignment from class UPtr which has operator-> defined resulting
   * in a value whcih can be const_cast'ed to T *.
   */
  template <typename UPtr>
  TransientRCPtr & assignConst(const UPtr & u) {
    ptr = const_cast<pointer>(PtrTraits<UPtr>::barePointer(u));
    return *this;
  }

  /**
   * Test for equality of the underlying pointers.
   */
  bool operator==(const TransientRCPtr & p) const { return ptr == p.ptr; }

  /**
   * Test for inequality of the underlying pointers.
   */
  bool operator!=(const TransientRCPtr & p) const { return ptr != p.ptr; }

  /**
   * Test for equality of the underlying pointers.
   */
  bool operator==(const_pointer p) const { return ptr == p; }

  /**
   * Test for inequality of the underlying pointers.
   */
  bool operator!=(const_pointer p) const { return ptr != p; }

  /**
   * Test for equality of the underlying pointers.
   */
  template <typename UPtr>
  bool operator==(const UPtr & u) const { return ptr == PtrTraits<UPtr>::barePointer(u); }

  /**
   * Test for inequality of the underlying pointers.
   */
  template <typename UPtr>
  bool operator!=(const UPtr & u) const { return ptr != PtrTraits<UPtr>::barePointer(u); }

  /**
   * Test for ordering.
   */
  bool operator<(const TransientRCPtr & p) const {
    return ( ptr && p.ptr && ptr->uniqueId != p.ptr->uniqueId ) ?
      ptr->uniqueId < p.ptr->uniqueId : ptr < p.ptr;
  }

  /**
   * Test for ordering.
   */
  bool operator<(const_pointer p) const {
    return ( ptr && p && ptr->uniqueId != p->uniqueId ) ?
      ptr->uniqueId < p->uniqueId : ptr < p;
  }

  /**
   * Returns true if the underlying pointer is null.
   */
  bool operator!() const { return !ptr; }

  /**
   * Returns the underlying pointer.
   */
  operator T * () const { return ptr; }

  /**
   * Member access.
   */
  pointer operator->() const { return ptr; }

  /**
   * Dereferencing.
   */
  reference operator*() const { return *ptr; }

private:

  /**
   * The actual pointer.
   */
  pointer ptr;

};

/**
 * TransientConstRCPtr is a simple wrapper around a bare const pointer
 * which can be assigned to and from an RCPtr and ConstRCPtr without
 * problem.
 *
 * @see RCPtr
 * @see ConstRCPtr
 */
template <typename T>
class TransientConstRCPtr {

public:

  /** Template argument typedef. */
  typedef void iterator_category;
  /** Template argument typedef. */
  typedef void difference_type;
  /** Template argument typedef. */
  typedef T * pointer;
  /** Template argument typedef. */
  typedef const T * const_pointer;
  /** Template argument typedef. */
  typedef T & reference;
  /** Template argument typedef. */
  typedef const T & const_reference;
  /** Template argument typedef. */
  typedef T value_type;

public:

  /**
   * Constructor for null pointer.
   */
  TransientConstRCPtr() : ptr(nullptr) {}

  /**
   * Constructor for nullptr.
   */
  TransientConstRCPtr( decltype(nullptr) ) : ptr(nullptr) {}

  /**
   * Copy constructor.
   */
  TransientConstRCPtr(const TransientConstRCPtr & p) : ptr(p.ptr) {}

  /**
   * Copyconstructor for class UPtr which has operator-> defined
   * resulting in a value implicitly convertible to const T *.
   */
  template <typename UPtr>
  TransientConstRCPtr(const UPtr & u) : ptr(PtrTraits<UPtr>::barePointer(u)) {}

  /**
   * Construction from real pointer.
   */
  explicit TransientConstRCPtr(const_pointer p) : ptr(p) {}

  /**
   * Destructor. Does not delete the object pointed to.
   */
  ~TransientConstRCPtr() {}

  /**
   * Assignment.
   */
  TransientConstRCPtr & operator=(const TransientConstRCPtr & p) {
    ptr = p.ptr;
    return *this;
  }

  /**
   * Assignment from class UPtr which has operator-> defined resulting
   * in a value implicitly convertible to const T *.
   */
  template <typename UPtr>
  TransientConstRCPtr & operator=(const UPtr & u) {
    ptr = PtrTraits<UPtr>::barePointer(u);
    return *this;
  }

  /**
   * Assignment from class UPtr which has operator-> defined resulting
   * in a value which can be cast dynamically to const T *.
   */
  template <typename UPtr>
  TransientConstRCPtr & assignDynamic(const UPtr & u) {
    ptr = dynamic_cast<const_pointer>(PtrTraits<UPtr>::barePointer(u));
    return *this;
  }

  /**
   * Test for equality of the underlying pointers.
   */
  bool operator==(const TransientConstRCPtr & p) const { return ptr == p.ptr; }

  /**
   * Test for inequality of the underlying pointers.
   */
  bool operator!=(const TransientConstRCPtr & p) const { return ptr != p.ptr; }

  /**
   * Test for equality of the underlying pointers.
   */
  bool operator==(const_pointer p) const { return ptr == p; }


  /**
   * Test for inequality of the underlying pointers.
   */
  bool operator!=(const_pointer p) const { return ptr != p; }

  /**
   * Test for equality of the underlying pointers.
   */
  template <typename UPtr>
  bool operator==(const UPtr & u) const { return ptr == PtrTraits<UPtr>::barePointer(u); }

  /**
   * Test for inequality of the underlying pointers.
   */
  template <typename UPtr>
  bool operator!=(const UPtr & u) const { return ptr != PtrTraits<UPtr>::barePointer(u); }

  /**
   * Test for ordering.
   */
  bool operator<(const TransientConstRCPtr & p) const {
    return ( ptr && p.ptr && ptr->uniqueId != p.ptr->uniqueId ) ?
      ptr->uniqueId < p.ptr->uniqueId : ptr < p.ptr;
  }

  /**
   * Test for ordering.
   */
  bool operator<(const_pointer p) const {
    return ( ptr && p && ptr->uniqueId != p->uniqueId ) ?
      ptr->uniqueId < p->uniqueId : ptr < p;
  }

  /**
   * Returns true if the underlying pointer is null.
   */
  bool operator!() const { return !ptr; }

  /**
   * Returns (not) the underlying pointer.
   */
  operator const T * () const { return ptr; }

  /**
   * Member access.
   */
  const_pointer operator->() const { return ptr; }

  /**
   * Dereferencing.
   */
  const_reference operator*() const { return *ptr; }

private:

  /**
   * The actual pointer.
   */
  const_pointer ptr;

};

/**
 * Specialization of the PtrTraits class for RCPtr.
 */
template <typename T>
struct PtrTraits< RCPtr<T> >: public PtrTraitsType {

  /** Template argument typedef. */
  typedef typename RCPtr<T>::value_type value_type;
  /** Template argument typedef. */
  typedef typename RCPtr<T>::reference reference;
  /** Template argument typedef. */
  typedef typename RCPtr<T>::const_reference const_reference;
  /** Template argument typedef. */
  typedef RCPtr<T> pointer;
  /** Template argument typedef. */
  typedef ConstRCPtr<T> const_pointer;
  /** Template argument typedef. */
  typedef TransientRCPtr<T> transient_pointer;
  /** Template argument typedef. */
  typedef TransientConstRCPtr<T> transient_const_pointer;

  /**
   * Return the bare pointer of the given pointer object.
   */
  static T * barePointer(const RCPtr<T> & p) { return p.operator->(); }

  /**
   * Create an object and return a pointer to it.
   */
  static pointer create() { return RCPtr<T>::Create(); }

  /**
   * Create an copy of an object and return a pointer to it.
   */
  static pointer create(const_reference t) { return RCPtr<T>::Create(t); }

  /**
   * Destroy the object pointed to.
   */
  static void destroy(pointer) {}

  /**
   * Cast dynamically.
   */
  template <typename UPtr>
  static pointer DynamicCast(const UPtr & u) {
    pointer t;
    t.assignDynamic(u);
    return t;
  }

  /**
   * Cast away constness.
   */
  template <typename UPtr>
  static pointer ConstCast(const UPtr & u) {
    pointer t;
    t.assignConst(u);
    return t;
  }

  /**
   * Cast from a basic pointer.
   */
  static pointer PtrCast(T * t) {
    return pointer(t);
  }

  /**
   * RCPtr is reference counted.
   */
  static const bool reference_counted = true;

}; 

/**
 * Specialization of the PtrTraits class for ConstRCPtr.
 */
template <typename T>
struct PtrTraits< ConstRCPtr<T> >: public PtrTraitsType {

  /** Template argument typedef. */
  typedef typename ConstRCPtr<T>::value_type value_type;
  /** Template argument typedef. */
  typedef typename ConstRCPtr<T>::reference reference;
  /** Template argument typedef. */
  typedef typename ConstRCPtr<T>::const_reference const_reference;
  /** Template argument typedef. */
  typedef RCPtr<T> pointer;
  /** Template argument typedef. */
  typedef ConstRCPtr<T> const_pointer;
  /** Template argument typedef. */
  typedef TransientRCPtr<T> transient_pointer;
  /** Template argument typedef. */
  typedef TransientConstRCPtr<T> transient_const_pointer;

  /**
   * Return the bare pointer of the given pointer object.
   */
  static const T * barePointer(const ConstRCPtr<T> & p) {
    return p.operator->();
  }

  /**
   * Create an object and return a pointer to it.
   */
  static const_pointer create() {
    return RCPtr<T>::Create();
  }

  /**
   * Create an copy of an object and return a pointer to it.
   */
  static const_pointer create(const_reference t) {
    return RCPtr<T>::Create(t);
  }

  /**
   * Destroy the object pointed to.
   */
  static void destroy(const_pointer) {}

  /**
   * Cast dynamically.
   */
  template <typename UPtr>
  static const_pointer DynamicCast(const UPtr & u) {
    const_pointer t;
    t.assignDynamic(u);
    return t;
  }

  /**
   * Cast away constness.
   */
  template <typename UPtr>
  static const_pointer ConstCast(const UPtr & u) {
    const_pointer t;
    t.assignDynamic(u);
    return t;
  }

  /**
   * Cast from a basic pointer.
   */
  static const_pointer PtrCast(const T * t) {
    return const_pointer(t);
  }

  /**
   * ConstRCPtr is reference counted.
   */
  static const bool reference_counted = true;

}; 

/**
 * Specialization of the PtrTraits class for TransientRCPtr.
 */
template <typename T>
struct PtrTraits< TransientRCPtr<T> >: public PtrTraitsType {

  /** Template argument typedef. */
  typedef typename TransientRCPtr<T>::value_type value_type;
  /** Template argument typedef. */
  typedef typename TransientRCPtr<T>::reference reference;
  /** Template argument typedef. */
  typedef typename TransientRCPtr<T>::const_reference const_reference;
  /** Template argument typedef. */
  typedef RCPtr<T> pointer;
  /** Template argument typedef. */
  typedef ConstRCPtr<T> const_pointer;
  /** Template argument typedef. */
  typedef TransientRCPtr<T> transient_pointer;
  /** Template argument typedef. */
  typedef TransientConstRCPtr<T> transient_const_pointer;

  /**
   * Return the bare pointer of the given pointer object.
   */
  static T * barePointer(const TransientRCPtr<T> & p) {
    return p.operator->();
  }

  /**
   * Destroy the object pointed to.
   */
  static void destroy(transient_pointer) {}

  /**
   * Cast dynamically.
   */
  template <typename UPtr>
  static transient_pointer DynamicCast(const UPtr & u) {
    transient_pointer t;
    t.assignDynamic(u);
    return t;
  }

  /**
   * Cast away constness.
   */
  static transient_pointer ConstCast(transient_const_pointer c) {
    transient_pointer t;
    t.assignConst(c);
    return t;
  }

  /**
   * Cast from a basic pointer.
   */
   static transient_pointer PtrCast(T * t) {
    return transient_pointer(t);
  }

  /**
   * TransientRCPtr is not reference counted.
   */
  static const bool reference_counted = false;

}; 

/**
 * Specialization of the PtrTraits class for TransientConstRCPtr.
 */
template <typename T>
struct PtrTraits< TransientConstRCPtr<T> >: public PtrTraitsType {

  /** Template argument typedef. */
  typedef typename TransientConstRCPtr<T>::value_type value_type;
  /** Template argument typedef. */
  typedef typename TransientConstRCPtr<T>::reference reference;
  /** Template argument typedef. */
  typedef typename TransientConstRCPtr<T>::const_reference const_reference;
  /** Template argument typedef. */
  typedef RCPtr<T> pointer;
  /** Template argument typedef. */
  typedef ConstRCPtr<T> const_pointer;
  /** Template argument typedef. */
  typedef TransientRCPtr<T> transient_pointer;
  /** Template argument typedef. */
  typedef TransientConstRCPtr<T> transient_const_pointer;

  /**
   * Return the bare pointer of the given pointer object.
   */
  static const T * barePointer(const TransientConstRCPtr<T> & p) {
    return p.operator->();
  }

  /**
   * Destroy the object pointed to.
   */
  static void destroy(transient_const_pointer) {}

  /**
   * Cast dynamically.
   */
  template <typename UPtr>
  static transient_const_pointer DynamicCast(const UPtr & u) {
    transient_const_pointer t;
    t.assignDynamic(u);
    return t;
  }

  /**
   * Cast away constness.
   */
  template <typename UPtr>
  static transient_const_pointer ConstCast(const UPtr & u) {
    transient_const_pointer t;
    t.assignConst(u);
    return t;
  }

  /**
   * Cast from a basic pointer.
   */
  static transient_const_pointer PtrCast(const T * t) {
    return transient_const_pointer(t);
  }

  /**
   * TransientConstRCPtr is not reference counted.
   */
  static const bool reference_counted = false;

}; 

}
}

namespace std {

/**
 * Specialization of std::swap to avoid unnecessary (in/de)crements of
 * the reference count.
 */
template <typename T>
inline void swap(ThePEG::Pointer::RCPtr<T> & t1, 
		 ThePEG::Pointer::RCPtr<T> & t2) {
  t1.swap(t2);
}

/**
 * Specialization of std::swap to avoid unnecessary (in/de)crements of
 * the reference count.
 */
template <typename T>
inline void swap(ThePEG::Pointer::ConstRCPtr<T> & t1,
		 ThePEG::Pointer::ConstRCPtr<T> & t2) {
  t1.swap(t2);
}

}

#endif /* ThePEG_RCPtr_H */
