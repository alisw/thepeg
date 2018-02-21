// -*- C++ -*-
//
// Rebinder.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Rebinder_H
#define ThePEG_Rebinder_H

#include "ThePEG/Config/ThePEG.h"
#include "Rebinder.fh"
#include <stdexcept>

namespace ThePEG {

/**
 * Rebinder is a class associating pairs of pointers to objects. It is
 * typically used when cloning a set of objects which have pointers to
 * eachother. The Rebinder is then set up so that a cloned object can
 * easily be retrieved given a pointer to the original one. The cloned
 * objects can then use the Rebinder to change their pointers so that
 * they henceforth points to the cloned copies.
 */
template <typename T>
class Rebinder {

public:
  
  /** Default constructor. */
  Rebinder() {}

public:

  ThePEG_DECLARE_TEMPLATE_POINTERS(T,TPtr);

  /** Map associating pairs of objects. */
  typedef std::map<cTPtr,TPtr> MapType;

  /** The iterator of the underlying map. */
  typedef typename MapType::const_iterator const_iterator;

public:

  /**
   * Return a pointer to the object associated with the argument.
   */
  TPtr & operator[](tcTPtr t) { return theMap[t]; }

  /**
   * Return a pointer to the object associated with the argument. If
   * no corresponding object is found a null pointer given by R() is
   * returned.
   * @param r a pointer to an object of a type which is derived from T.
   */
  template <typename R>
  R translate(const R & r) const {
    const_iterator it = theMap.find(r);
    return it == theMap.end()? R(): dynamic_ptr_cast<R>(it->second);
  }

  /**
   * Insert pointers to objects into the output iterator, pointers to
   * objects corresponding to the ones given by the range of input
   * iterators. If a given object in the input iterator range does not
   * exists, a null pointer will be inserted in the output iterator.
   */
  template <typename OutputIterator, typename InputIterator>
  void translate(OutputIterator r,
		 InputIterator first, InputIterator last) const {
    while ( first != last ) *r++ = translate(*first++);
  }

  /**
   * Return a pointer to the object associated with the argument. If
   * no corresponding object is found an exception is thrown.
   * @param r a pointer to an object of a type which is derived from T.
   */
  template <typename R>
  R alwaysTranslate(const R & r) const {
    R ret;
    if ( !r ) return ret;
    const_iterator it = theMap.find(r);
    ret = (it == theMap.end()? R(): dynamic_ptr_cast<R>(it->second));
    if ( !ret ) throw std::runtime_error("Rebinder exception");
    return ret;
  }

  /**
   * Insert pointers to objects into the output iterator, pointers to
   * objects corresponding to the ones given by the range of input
   * iterators. If a given object in the input iterator range does not
   * exists, an exception will be thrown.
   */
  template <typename OutputIterator, typename InputIterator>
  void alwaysTranslate(OutputIterator r, InputIterator first, InputIterator last)
    const {
    while ( first != last ) *r++ = alwaysTranslate(*first++);
  }

  /**
   * Acces the underlying map representation.
   */
  const MapType & map() const { return theMap; }

private:


  /**
   * The underlying map representation.
   */
  MapType theMap;

private:

  /** The copy constructor is private and not implemented */
  Rebinder(const Rebinder &);

  /** The assignment operator is private and not implemented */
  Rebinder & operator=(const Rebinder &);

};


}

#endif /* ThePEG_Rebinder_H */
