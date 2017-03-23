// -*- C++ -*-
//
// Triplet.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Triplet_H
#define ThePEG_Triplet_H

#include "ThePEG/Config/ThePEG.h"

namespace ThePEG {

/**
 * The Triplet class represents a general triplet of objects
 * completely analogous to <code>std::pair</code>.
 */
template <typename T1, typename T2, typename T3>
struct Triplet {

  /** The type of the first member. */
  typedef T1 first_type;
  /** The type of the second member. */
  typedef T2 second_type;
  /** The type of the third member. */
  typedef T3 third_type;

  /** The first member. */
  T1 first;
  /** The second member. */
  T2 second;
  /** The third member. */
  T3 third;

  /** Default construcotr. */
  Triplet() : first(T1()), second(T2()), third(T3()) {}

  /** Constructor specifying the three members. */
  Triplet(const T1 & t1, const T2 & t2, const T3 & t3)
    : first(t1), second(t2), third(t3) {}

  /** Copy constructor. */
  Triplet(const Triplet<T1,T2,T3> & t)
    : first(t.first), second(t.second), third(t.third) {}

  /** Copy constructor from other Triplet type. */
  template <typename U1, typename U2, typename U3>
  Triplet(const Triplet<U1,U2,U3> & u)
    : first(u.first), second(u.second), third(u.third) {}

  /** Test for equality. */
  bool operator==(const Triplet<T1,T2,T3> & t) const {
    return first == t.first && second == t.second && third == t.third;
  }

  /** Test for ordering.
   * @return <code>first < t.first || ( * !(t.first < first) && (
   * second < t.second || ( !(t.second < second) && third < t.third
   * )))</code> */
  bool operator<(const Triplet<T1,T2,T3> & t) const {
    return first < t.first ||
      ( !(t.first < first) && 
	( second < t.second ||
	  ( !(t.second < second) && third < t.third )));
  }
};

/** Helper function returning a Triplet with template parameters
 *  determined by the arguments. */
template <typename T1, typename T2, typename T3>
inline Triplet<T1,T2,T3>
makeTriplet (const T1 & t1, const T2 & t2, const T3 & t3)
{
  return Triplet<T1,T2,T3>(t1, t2, t3);
}

/** Output a Triplet to a stream. */
template <typename OStream, typename T1, typename T2, typename T3>
OStream & operator<<(OStream & os, const Triplet<T1,T2,T3> & t) {
  return os << t.first << t.second << t.third;
}

/** Input a Triplet from a stream. */
template <typename IStream, typename T1, typename T2, typename T3>
IStream & operator>>(IStream & is, Triplet<T1,T2,T3> & t) {
  return is >> t.first >> t.second >> t.third;
}

}

#endif /* ThePEG_Triplet_H */
