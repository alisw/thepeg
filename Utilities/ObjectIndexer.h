// -*- C++ -*-
//
// ObjectIndexer.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_ObjectIndexer_H
#define THEPEG_ObjectIndexer_H
// This is the declaration of the ObjectIndexer class.

#include "ThePEG/Config/ThePEG.h"
#include <limits>

namespace ThePEG {

/**
 * This is a templated class which dynamically associates (reference
 * counted) objects to integer indices. By default, all indices will
 * be non-negative, but explicit usage of negative indices is
 * allowed as long as they do not include NoIndex.
 */
template <typename IntT, typename ObjT, IntT NoIndex = static_cast<IntT>(-1)>
class ObjectIndexer {

public:

  ThePEG_DECLARE_TEMPLATE_POINTERS(ObjT,TPtr);

  /** Map of objects to indices */
  typedef map<IntT,tTPtr> IndexObjectMap;

  /** Map of indices to objects. */
  typedef map<TPtr,IntT> ObjectIndexMap;

public:

  /**
   * Empty constructor.
   */
  ObjectIndexer(): next(0) {}

  /**
   * Return the index for the given object. If the object is not known,
   * a new index will be created.
   */
  IntT operator()(tTPtr o) {
    typename ObjectIndexMap::iterator it = objectIndex.find(o);
    if ( it == objectIndex.end() ) {
      IntT i = next++;
      objectIndex[o] = i;
      indexObject[i] = o;
      return i;
    } else
      return it->second;
  }

  /**
   * Return the index for the given object. If the object is not known,
   * NoIndex will be returned.
   */
  IntT operator()(tTPtr o) const {
    return find(o);
  }

  /**
   * Return the index for the given object. If the object is not known,
   * NoIndex will be returned.
   */
  IntT find(tTPtr o) const {
    typename ObjectIndexMap::const_iterator it = objectIndex.find(o);
    return it == objectIndex.end()? NoIndex: it->second;
  }

  /**
   * Return the object for the given index. If the index is not known,
   * a new object will be (default) created.
   */
  tTPtr operator()(IntT i) {
    if ( i == NoIndex ) return tTPtr();
    typename IndexObjectMap::iterator it = indexObject.find(i);
    if ( it == indexObject.end() ) {
      TPtr o = new_ptr<ObjT>();
      objectIndex[o] = i;
      indexObject[i] = o;
      next = max(next, i + 1);
      return o;
    } 
    else
      return it->second;
  }

  /**
   * Return the object for the given index. If the index is not known,
   * a null pointer will be returned.
   */
  tTPtr operator()(IntT i) const {
    return find(i);
  }

  /**
   * Return the object for the given index. If the index is not known,
   * a null pointer will be returned.
   */
  tTPtr find(IntT i) const {
    typename IndexObjectMap::const_iterator it = indexObject.find(i);
    return it == indexObject.end()? tTPtr(): it->second;
  }

  /**
   * Associate the given object with the given index. Possible other
   * associations involving the index or the object is removed. If the
   * given index is NoIndex, this function does nothing.
   */
  void operator()(IntT i, tTPtr o) {
    if ( i == NoIndex ) return;
    typename IndexObjectMap::iterator iit = indexObject.find(i);
    if ( iit != indexObject.end() ) objectIndex.erase(iit->second);
    typename ObjectIndexMap::iterator oit = objectIndex.find(o);
    if ( oit != objectIndex.end() ) indexObject.erase(oit->second);
    objectIndex[o] = i;
    indexObject[i] = o;
    next = max(next, i + 1);
  }

  /**
   * Return true if the given object is known.
   */
  bool included(tTPtr o) const {
    return objectIndex.find(o) != objectIndex.end();
  }

  /**
   * Return true if the given index is known.
   */
  bool included(IntT i) const {
    return indexObject.find(i) != indexObject.end();
  }

  /**
   * Remove all associations.
   */
  void clear() {
    indexObject.clear();
    objectIndex.clear();
  }

  /**
   * Return true if no associations has been made.
   */
  bool empty() const {
    return indexObject.empty() && objectIndex.empty();
  }

private:

  /**
   * All known objects keyed by their indices.
   */
  IndexObjectMap indexObject;

  /**
   * All known indices keyed by the corresponding objects.
   */
  ObjectIndexMap objectIndex;

  /**
   * The next index to be used.
   */
  IntT next;

private:

  /**
   * Private and non-existent assignment operator.
   */
  ObjectIndexer & operator=(const ObjectIndexer &) = delete;

};

}

#endif /* THEPEG_ObjectIndexer_H */
