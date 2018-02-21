// -*- C++ -*-
//
// ReferenceCounted.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ReferenceCounted_H
#define ThePEG_ReferenceCounted_H
// This is the declaration of the ReferenceCounted class.


#include "RCPtr.fh"
#include "ThePEG/Persistency/PersistentIStream.fh"

namespace ThePEG {
namespace Pointer {

/**
 * ReferenceCounted must be the (virtual) base class of all
 * classes which may be pointed to by the RCPtr smart
 * pointer class. It keeps track of all RCPtr and
 * ConstRCPtr pointers which are currently pointing to an
 * object.
 *
 * @see RCPtr
 * @see ConstRCPtr
 */
class ReferenceCounted {

  /** The RCPtrBase class needs to acces the private parts of ReferenceCounted. */
  friend class RCPtrBase;
  friend class ThePEG::PersistentIStream;

public:

  /**
   * The integer type used for counting.
   */
  typedef unsigned int CounterType;

protected:

  /** @name Standard constructors and assignment. */
  //@{
  /**
   * Default constructor.
   */
  ReferenceCounted() 
    : uniqueId(++objectCounter), 
      theReferenceCounter(CounterType(1)) {}

  /**
   * Copy-constructor.
   */
  ReferenceCounted(const ReferenceCounted &)
    : uniqueId(++objectCounter), 
      theReferenceCounter(CounterType(1)) {}

  /**
   * Assignment.
   */
  ReferenceCounted & operator=(const ReferenceCounted &)
  {
    return *this;
  }
    
  //@}

public:

  /**
   * Return the reference count.
   */
  CounterType referenceCount() const 
  { 
    return theReferenceCounter; 
  }

private:

  /**
   * Increment the reference count.
   */
  void incrementReferenceCount() const 
  { 
    ++theReferenceCounter; 
  }

  /**
   * Decrement with the reference count.
   */
  bool decrementReferenceCount() const 
  {
    return !--theReferenceCounter;
  }

public:

  /**
   * The unique ID. Can be used as sorting criterion, e.g. for set<> and map<>.
   */
  const unsigned long uniqueId;

private:

  /** 
   * A counter for issuing unique IDs. It will overflow back to 0 eventually,
   * but it is very unlikely that two identical IDs show up in the same event.
   */
  static unsigned long objectCounter;

  /**
   * The reference count.
   */
  mutable CounterType theReferenceCounter;

};


}
}

#endif /* ThePEG_ReferenceCounted_H */

