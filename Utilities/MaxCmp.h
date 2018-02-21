// -*- C++ -*-
//
// MaxCmp.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_MaxCmp_H
#define THEPEG_MaxCmp_H
//
// This is the declaration of the MaxCmp class.
//

#include <functional>

namespace ThePEG {

/**
 * MaxCmp is a helper class to be used in a loop where one would like
 * to keep track of the largest value so far of a certain
 * expression. The class simply checks if the given value to the
 * operator() is the largest so far (in which case true is returned,
 * and the value is saved together with the optional index
 * argument. The Cmp template argument is by default greater<T>, but
 * can be set to any comparison class to change the meaning of
 * maximum: MaxCmp<double, int, less<double> > will keep track of the
 * smallest value.
 */
template <typename T = double, typename Indx = int, typename Cmp = std::greater<T> >
class MaxCmp {

public:

  /**
   * The default constructor.
   */
  MaxCmp() : init(false), max(T()), indx(Indx()) {}

  /**
   * Constructor specifying an initial maximum value, \a t.
   */
  MaxCmp(const T & t, Indx in = Indx()) : init(true), max(t), indx(in) {}

public:

  /**
   * If \a t is the largest value seen so far return true. Otherwise
   * return false. \a i is an optional index for the value \a t.
   */
  bool operator()(const T & t, Indx i = Indx())
  {
    if ( !init || cmp(t, max) ) {
      max = t;
      init = true;
      indx = i;
      return true;
    }
    return false;
  }

  /**
   * Return the largest value so far.
   */
  operator const T & () const { return value(); }

  /**
   * Return the largest value so far.
   */
  const T & value() const { return max; }

  /**
   * Return the index of the largest object seen so far.
   */
  Indx index() const {
    return indx;
  }

  /**
   * Return true if no index has been chosen.
   */
  bool operator!() const {
    return !init;
  }

private:

  /**
   * True if a first value has been given;
   */
  bool init;

  /**
   * The largest value seen so far.
   */
  T max;

  /**
   * The index for the largest value seen so far.
   */
  Indx indx;

  /**
   * The comparison object to be used.
   */
  Cmp cmp;

};

/**
 * Special calss for Minimum comparisons.
 */
template <typename T, typename Indx = int>
class MinCmp: public MaxCmp<T, Indx, less<T> > {

public:

  /**
   * Constructors are not inherited.
   */
  MinCmp() {}
  
  /**
   * Constructors are not inherited.
   */
  MinCmp(const T & t, Indx in = Indx()) : MaxCmp<T, Indx, less<T> >(t, in) {}
  
};

}

#endif /* THEPEG_MaxCmp_H */
