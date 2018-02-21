// -*- C++ -*-
//
// Interval.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Interval_H
#define ThePEG_Interval_H
// This is the declaration of the Interval class.

#include <utility>
#include <vector>
#include "Interval.fh"
#include "ThePEG/Utilities/UnitIO.h"

namespace ThePEG {

template <typename T, typename CMP>
/**
 * An <code>Interval</code> object is used to represent an interval
 * <code>[ lower(), upper() )</code> where the ordering is defined by
 * the <code>bool operator()(const T &, const T &) const</code> member
 * of the CMP class (by defaut less<T>).
 */
class Interval {

public:

  /**
   * Construct an empty interval.
   */
  Interval() : theLimits(pair<T,T>()) {}

  /**
   * Construct interval [dn,up).
   */
  Interval(T dn, T up) : theLimits(pair<T,T>(dn, up)) {}

  /**
   * Test for equality.
   */
  bool operator==(const Interval & i) const {
    return lower() == i.lower() && upper() == i.upper();
  }

  /**
   * Test for ordering.
   * @return true if <code>lower() < i.lower()</code> or <code>lower()
   * == i.lower()</code> and <code>upper() < i.upper()</code>.
   */
  bool operator<(const Interval & i) const {
    return cmp(lower(), i.lower()) ||
      ( lower() == i.lower() && cmp(upper(), i.upper()) );
  }

  /**
   * Check consistency ie. that lower() < upper().
   */
  bool check() const { return cmp(lower(), upper()); }

  /**
   * Returns true if x is within the interval.
   */
  bool operator()(T x) const { return includes(x); }

  /**
   * Returns true if x is within the interval.
   */
  bool includes(T x) const { return !cmp(x, lower()) && cmp(x, upper()); }

  /**
   * Returns true if the whole of i is within the interval.
   */
  bool includes(const Interval<T,CMP> & i) const {
    return includes(i.lower()) && !cmp(upper(), i.upper());
  }

  /**
   * If x is in the interval return the interval [x,upper()) and
   * change this interval to [lower(),x). If x is not within the
   * interval, return [0,0) and leave this interval as it is.
   */
  Interval<T,CMP> chopUpper(T x) {
    Interval<T,CMP> r;
    if ( includes(x) ) {
      r.lower(x);
      r.upper(upper());
      upper(x);
    }
    return r;
  }

  /**
   * If x is in the interval return the interval [lower(),x) and
   * change this interval to [x,upper()). If x is not within the
   * interval, return [0,0) and leave this interval as it is.
   */
  Interval<T,CMP> chopLower(T x) {
    Interval<T,CMP> r;
    if ( includes(x) ) {
      r.lower(lower());
      r.upper(x);
      lower(x);
    }
    return r;
  }

  /**
   * If this interval operlaps with i return the overlapping interval.
   */
  Interval<T,CMP> overlap(const Interval & i) const {
    Interval<T,CMP> res;
    if ( operator==(i) ) res = i;
    if ( includes(i.upper()) || includes(i.lower()) )
      res = Interval<T,CMP>(max(lower(),i.lower()), min(upper(), i.upper()));
    return res;
  }

  /**
   * If this interval operlaps with i return the union of the two
   * intervals.
   */
  Interval<T,CMP> sum(const Interval & i) const {
    Interval<T,CMP> res;
    if ( operator==(i) ) res = i;
    if ( includes(i.upper()) || includes(i.lower()) )
      res = Interval<T,CMP>(min(lower(),i.lower()), max(upper(), i.upper()));
    return res;
  }

  /**
   * Return the upper limit of the interval.
   */
  T upper() const { return theLimits.second; }

  /**
   * Return the lower limit of the interval.
   */
  T lower() const { return theLimits.first; }

  /**
   * Set the upper limit of the interval.
   */
  void upper(T up) { theLimits.second = up; }

  /**
   * Set the lower limit of the interval.
   */
  void lower(T dn) { theLimits.first = dn; }

  /**
   * Check if any of the values in the iterator given range is
   * included in this interval.
   */
  template <typename Iterator>
  bool check(Iterator first, Iterator last);

  /**
   * Check if all of the values in the given iterator range is
   * included in this interval.
   */
  template <typename Iterator>
  bool checkAll(Iterator first, Iterator last);

  /**
   * If x is in the given interval, split the given interval in two,
   * otherwise return an empty list.
   */
  std::vector< Interval<T,CMP> > split(Interval<T,CMP>, T x);
  
  /**
   * For each value in the given range is in the given interval, split
   * the interval in two, otherwise return an empty list.
   */
  template<typename Iterator>
  std::vector< Interval<T,CMP> > split(Interval<T,CMP>,
				       Iterator first, Iterator last);

private:

  /** The lower and upper limit of this interval */
  std::pair<T,T> theLimits;

  /** The object used for comparisons. */
  static CMP cmp;

};

/** An interval of doubles. */
typedef Interval<double> DInterval;

/** Helper function to create an interval of a type determined by the
 *  parameters. */
template <typename T, typename CMP>
inline Interval<T,CMP> makeInterval(T dn, T up) { return Interval<T,CMP>(dn, up); }

/** Ouptut an interval to a stream. */
template <typename OStream, typename T, typename CMP>
inline OStream & operator<<(OStream & os, const Interval<T,CMP> & i) {
  os << i.lower() << i.upper();
  return os;
}

/** Input an interval from a stream. */
template <typename IStream, typename T, typename CMP>
inline IStream & operator>>(IStream & is, Interval<T,CMP> & i) {
  T up, dn;
  is >> dn >> up;
  i.lower(dn);
  i.upper(up);
  return is;
}

/** Output an interval of a diminsionful type to a stream using the
 *  given unit.
 * @param os the output stream.
 * @param i the interval.
 * @param u the unit. */
template <typename OStream, typename T, typename CMP, typename UT>
void ounitstream(OStream & os, const Interval<T,CMP> & i, UT & u) {
  os << ounit(i.lower(), u) << ounit(i.upper(), u);
}

/** Input an interval of a diminsionful type from a stream using the
 *  given unit.
 * @param is the input stream.
 * @param i the interval.
 * @param u the unit. */
template <typename IStream, typename T, typename CMP, typename UT>
void iunitstream(IStream & is, Interval<T,CMP> & i, UT & u) {
  T low, upp;
  is >> iunit(low, u) >> iunit(upp, u);
  i = Interval<T,CMP>(low, upp);
}

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "Interval.tcc"
#endif

#endif /* ThePEG_Interval_H */
