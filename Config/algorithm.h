// -*- C++ -*-
//
// algorithm.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_algorithm_H
#define ThePEG_algorithm_H

/** \file
 * This file implements a number of interfaces to <code>std::</code>
 * algorithms, modified to take a whole container as argument rather
 * than a range of iterators, Also defines IteratorRange to
 * encapsulate a range of iterators and corresponding algorithms.
 */

#include "ThePEG/Config/ThePEG.h"
#include <algorithm>

namespace ThePEG {

/**
 * A pair of iterators to be used in specialized algorithms instead
 * of the standard first, last construction.
 */
template <typename Iterator>
struct IteratorRange: public std::pair<Iterator,Iterator> {

  /** The underlying representation. */
  typedef std::pair<Iterator,Iterator> BaseType;

  /** Default constructor. */
  IteratorRange() {}

  /** Copy constructor */
  IteratorRange(const IteratorRange & ir): BaseType(ir) {}

  /** Constructor taking the underlying pair representation as
      argument. */
  IteratorRange(const BaseType & ir): BaseType(ir) {}

};

/** Return an IteratorRange corresponding to the whole container. */
template <typename Container>
inline IteratorRange<typename Container::iterator>
range(Container & c) {
  return std::make_pair(c.begin(), c.end());
}

/** Return an IteratorRange of const iterators corresponding to the
 *  whole container. */
template <typename Container>
inline IteratorRange<typename Container::const_iterator>
range(const Container & c) {
  return std::make_pair(c.begin(), c.end());
}

/** Return an IteratorRange of reverse iterators corresponding to the
 *  whole container. */
template <typename Container>
inline IteratorRange<typename Container::reverse_iterator>
rrange(Container & c) {
  return std::make_pair(c.rbegin(), c.rend());
}

/** Return an IteratorRange of reverse const iterators corresponding
 *  to the whole container. */
template <typename Container>
inline IteratorRange<typename Container::const_reverse_iterator>
rrange(const Container & c) {
  return std::make_pair(c.rbegin(), c.rend());
}

/** The std::for_each function taking an IteratorRange as argument. */
template <typename Iterator, typename FNC>
inline FNC for_each(IteratorRange<Iterator> r, FNC f) {
  return std::for_each(r.first, r.second, f);
}

/** The std::find function taking an IteratorRange as argument. */
template <typename Iterator, typename T>
inline Iterator find(IteratorRange<Iterator> r, const T & t) {
  return std::find(r.first, r.second, t);
}

/** The std::find_if function taking an IteratorRange as argument. */
template <typename Iterator, typename Pred>
inline Iterator find_if(IteratorRange<Iterator> r, Pred p) {
  return std::find_if(r.first, r.second, p);
}

/** The std::replace function taking an IteratorRange as argument. */
template <typename Iterator, typename T>
inline void replace(IteratorRange<Iterator> r, const T & oval, const T & nval) {
  return std::replace(r.first, r.second, oval, nval);
}

/** The std::for_each function taking a whole container as argument. */
template <typename Cont, typename FNC>
inline FNC for_each(Cont & c, FNC f) {
  return std::for_each(c.begin(), c.end(), f);
}

/** The std::for_each function taking a whole const container as argument. */
template <typename Cont, typename FNC>
inline FNC for_each(const Cont & c, FNC f) {
  return std::for_each(c.begin(), c.end(), f);
}

/** The std::find function taking a whole container as argument. */
template <typename Cont, typename Type>
inline typename Cont::iterator find(Cont & c, const Type & t) {
  return find(range(c), t);
}

/** The std::find function taking a whole const container as argument. */
template <typename Cont, typename Type>
inline typename Cont::const_iterator find(const Cont & c, const Type & t) {
  return find(range(c), t);
}

/** The std::find_if function taking a whole container as argument. */
template <typename Cont, typename Pred>
inline typename Cont::iterator find_if(Cont & c, const Pred & p) {
  return find_if(range(c), p);
}

/** The std::find_if function taking a whole const container as argument. */
template <typename Cont, typename Pred>
inline typename Cont::const_iterator find_if(const Cont & c, const Pred & p) {
  return find_if(range(c), p);
}

/** The std::replace function taking a whole container as argument. */
template <typename Cont, typename T>
inline void replace(Cont & c, const T & oval, const T & nval) {
  return replace(range(c), oval, nval);
}

}

// #include "algorithm.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "algorithm.tcc"
#endif

#endif /* ThePEG_algorithm_H */
