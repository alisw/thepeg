// -*- C++ -*-
//
// DIterator.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_DIterator_H
#define ThePEG_DIterator_H
// This is the declaration of the DIterator class.

#include <iterator>
#include "ThePEG/Pointer/PtrTraits.h"

namespace ThePEG {

template <typename BaseIterator>
/**
 * DIterator is an iterator adaptor class. It can be used whenever one
 * has a container with pointers to facilitate member selection. The
 * only requirement is that the underlying pointer is pointing to a
 * valid object. Given e.g. a vector of pointers, <code>vector<A*>
 * pv</code> where the class <code>A</code> has a member function
 * <code>dosomething()</code>, it can be used as follows:<br>
 * <code>typedef DIterator<typename vector<A*>::iterator> It;</code><BR>
 * <code>for ( It i = dv.begin(), i != dv.end(), ++i )
 * i->dosomething();</code><BR>
 */
class DIterator {

public:

  /** Forward typedef from underlying iterator */
  typedef std::iterator_traits<BaseIterator>  	Traits;
  /** Forward typedef from underlying iterator */
  typedef typename Traits::iterator_category 	iterator_category;
  /** Forward typedef from underlying iterator */
  typedef typename Traits::difference_type 	difference_type;
  /** Forward typedef from underlying iterator */
  typedef typename Traits::value_type 		PtrType;
  /** Forward typedef from underlying iterator */
  typedef Pointer::PtrTraits<PtrType>         PtrTraits;
  /** Forward typedef from underlying iterator */
  typedef typename PtrTraits::value_type 	value_type;
  /** Forward typedef from underlying iterator */
  typedef typename PtrTraits::pointer          	pointer;
  /** Forward typedef from underlying iterator */
  typedef typename PtrTraits::reference         reference;

public:

  /**
   * Constructor from a normal iterator.
   */
  DIterator(const BaseIterator & in): i(in) {}

  /**
   * Copy constructor.
   */
  DIterator(const DIterator & pi): i(pi.i) {}
    
public:

  /**
   * Dereference the pointer referred to by the underlying iterator.
   */
  reference operator*() const { return **i; }

  /**
   * Select member from  the pointer referred to by the underlying iterator.
   */
  pointer operator->() const { return *i; }

  /**
   * Standard assignment operator.
   */
  DIterator & operator=(const DIterator & pi) { i = pi.i; return *this; }

  /**
   * Assignment from a a normal iterator.
   */
  DIterator & operator=(const BaseIterator & pi) { i = pi; return *this; }

  /** @name Increment and decrement operators. */
  //@{
  /** Pre increment the underlying iterator. */
  DIterator & operator++() { ++i; return *this; }
  /** Post increment the underlying iterator. */
  DIterator operator++(int) { DIterator tmp(*this); ++i; return tmp; }
  /** Pre decrement the underlying iterator. */
  DIterator & operator--() { --i; return *this; }
  /** Post decrement the underlying iterator. */
  DIterator operator--(int) { DIterator tmp(*this); --i; return tmp; }
  /** Jump forward n steps */
  DIterator & operator+=(int n) { i += n; return *this; }
  /** Get an iterator n steps forward. */
  DIterator operator+(int n) { return DIterator(i + n); }
  /** Jump backward n steps */
  DIterator & operator-=(int n) { i -= n; return *this; }
  /** Get an iterator n steps backward. */
  DIterator operator-(int n) { return DIterator(i - n); }
  //@}

  /**
   * Select a pointer with the given index and return a reference to
   * the object pointed to.
   */
  reference operator[](difference_type n) { return *(i[n]); }

  /**
   * Return the distance to the given iterator.
   */
  difference_type operator-(const DIterator & pi) { return i - pi.i; }

  /** @name Comparison operators. */
  //@{
  /** Test for equality. */
  bool operator==(const DIterator & pi) { return i == pi.i; }
  /** Test for inequality. */
  bool operator!=(const DIterator & pi) { return i != pi.i; }
  /** Test for less. */
  bool operator<(const DIterator & pi) { return i < pi.i; }
  /** Test for greater. */
  bool operator>(const DIterator & pi) { return i > pi.i; }
  /** Test for less or equal. */
  bool operator<=(const DIterator & pi) { return i <= pi.i; }
  /** Test for greater or equal. */
  bool operator>=(const DIterator & pi) { return i >= pi.i; }
  /** Test for equality. */
  bool operator==(const BaseIterator & bi) { return i == bi; }
  /** Test for inequality. */
  bool operator!=(const BaseIterator & bi) { return i != bi; }
  /** Test for less. */
  bool operator<(const BaseIterator & bi) { return i < bi; }
  /** Test for greater. */
  bool operator>(const BaseIterator & bi) { return i > bi; }
  /** Test for less or equal. */
  bool operator<=(const BaseIterator & bi) { return i <= bi; }
  /** Test for greater or equal. */
  bool operator>=(const BaseIterator & bi) { return i >= bi; }
  //@}

private:

  /**
   * The underlying standard iterator.
   */
  BaseIterator i;

private:

  /**
   * The default constructor should never be used.
   */
  DIterator() {}

};

}

#endif /* ThePEG_DIterator_H */
