// -*- C++ -*-
//
// ACDCTraits.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ACDCTraits_H
#define ACDCTraits_H

namespace ACDCGenerator {

/**
 * ACDCTraitsType is an empty non-polymorphic base class for all
 * traits classes in ACDCGenerator.
 */
struct ACDCTraitsType {};

/**
 * ACDCFncTraits defines the interface to functions to be sampled by
 * ACDCGen. It only defines one function which defines how the
 * functions are called. If the default implementation is not
 * suitable, ACDCFncTraits may be specialized for a function class
 * implementing a function with the same signature.
 */
template <typename FncPtr>
struct ACDCFncTraits: public ACDCTraitsType {

  /**
   * Call a function to be sampled by ACDCGen.
   * @return <code>(*f)(x)</code>.
   */
  static inline double value(const FncPtr & f, const DVector & x) {
    return (*f)(x);
  }

};

/**
 * ACDCRandomTraits defines the interface to random number generator
 * objects to be used by ACDCGen. If this default implementation is
 * not suitable, ACDCRandomTraits may be specialized for any class as
 * long as functions with the same signature are present.
 */
template <typename Rnd>
struct ACDCRandomTraits: public ACDCTraitsType {

  /**
   * Return a flat random number in the interval ]0,1[.
   */
  static inline double rnd(Rnd * r) { return r->flat(); }

  /**
   * Return a flat random number in the interval ]\a xl,\a xu[.
   */
  static inline double rnd(Rnd * r, double xl, double xu) {
    return xl + (xu - xl)*rnd(r);
  }

  /**
   * Generate a set of random numbers.
   * @param r the random generator.
   * @param l an input iterator giving the lower limit of the interval
   * of the first requested random number.
   * @param lend an input iterator marking the end of the range of
   * requested random numbers.
   * @param u an input iterator giving the upper limit of the interval
   * of the first requested random number.
   * @param res the ouput iterator used to output the random numbers.
   */
  template <typename InputIterator, typename OutputIterator>
  static inline void rnd(Rnd * r, InputIterator l, InputIterator lend,
			   InputIterator u, OutputIterator res) {
    for ( ; l != lend; ++l ) *res++ = *l + (*u++ - *l)*rnd(r);
  }

  /**
   * Generate \a D random numbers. The numbers are put into the
   * OutputIterator \a res.
   */
  template <typename OutputIterator>
  static inline void rnd(Rnd * r, int D, OutputIterator res) {
    for ( int d = 0; d < D; ++d ) *res++ = rnd(r);
  }

  /**
   * Return true with probability \a x.
   */
  static inline bool rndBool(Rnd * r, double x) { return rnd(r) < x; }

  /**
   * Return true with probability \a x(\a x + \a y).
   */
  static inline bool rndBool(Rnd * r, double x, double y) {
    return rndBool(r, x/(x + y)); }

  /**
   * Return a random integer in the interval [0,\a x[.
   */
  static inline long rndInt(Rnd * r, long x) {
    return long(rnd(r, 0.0, double(x)));
  }

};

}

#endif
