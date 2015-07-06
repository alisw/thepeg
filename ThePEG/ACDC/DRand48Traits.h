// -*- C++ -*-
//
// DRand48Traits.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef DRand48Traits_H
#define DRand48Traits_H

#include <cstdlib>
#include "ACDCTraits.h"

namespace ACDCGenerator {

/** @cond TRAITSPECIALIZATIONS */

/** Dummy struct to represent the standard drand48() random number generator. */
struct DRAND48 {};

/**
 * Specialization of ACDCRandomTraits for using the standard drand48()
 * random number generator.
 */
template <>
struct ACDCRandomTraits<DRAND48>: public ACDCTraitsType {

  /**
   * Return a flat random number in the interval ]0,1[.
   */
  static inline double rnd(DRAND48 * r) { return drand48(); }

  /**
   * Return a flat random number in the interval ]\a xl,\a xu[.
   */
  static inline double rnd(DRAND48 * r, double xl, double xu) {
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
  static inline void rnd(DRAND48 * r, InputIterator l, InputIterator lend,
			   InputIterator u, OutputIterator res) {
    for ( ; l != lend; ++l ) *res++ = *l + (*u++ - *l)*rnd(r);
  }

  /**
   * Generate \a D random numbers. The numbers are put into the
   * OutputIterator \a res.
   */
  template <typename OutputIterator>
  static inline void rnd(DRAND48 * r, int D, OutputIterator res) {
    for ( int d = 0; d < D; ++d ) *res++ = rnd(r);
  }

  /**
   * Return true with probability \a x.
   */
  static inline bool rndBool(DRAND48 * r, double x) { return rnd(r) < x; }

  /**
   * Return true with probability \a x(\a x + \a y).
   */
  static inline bool rndBool(DRAND48 * r, double x, double y) {
    return rndBool(r, x/(x + y)); }

  /**
   * Return a random integer in the interval [0,\a x[.
   */
  static inline long rndInt(DRAND48 * r, long x) {
    return long(rnd(r, 0.0, double(x))); }

};

/** @endcond */

}

#endif
