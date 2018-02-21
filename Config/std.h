// -*- C++ -*-
//
// std.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_std_H
#define ThePEG_std_H

/** \file

 * This file introduces a number of <code>std::</code> classes into
 * the ThePEG namespace. Also introduces some useful functions for
 * standard library classes.
 *
 * Do not make changes in this file. If you want to use alternatives
 * to the <code>std::</code> classes in ThePEG, edit a copy of this
 * file and include it in an alternative config file which can be
 * included in the main ThePEG.h config file using the macro
 * <code>ThePEG_ALTERNATE_CONFIG</code>.
 */

#include <map>
#include <set>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stack>
#include <utility>
#include <typeinfo>
#include <stdexcept>
#include <cmath>

namespace std {

/** @cond TRAITSPECIALIZATIONS */
/**
 * This specialization of the std::less class is needed in order to be
 * able use put pointers to type_info objects as keys in maps and
 * sets.
 */
template <>
struct less<const type_info *> :
    public binary_function<const type_info *, const type_info *, bool> 
{
  /**
   * This is the function called when comparing two pointers to
   * type_info.
   */
  bool operator()(const type_info * x, const type_info * y) const {
    return x->before(*y); }
};
/** @endcond */

}

namespace ThePEG {

using std::deque;
using std::stack;
using std::vector;
using std::multiset;
using std::set;
using std::map;
using std::list;
using std::multimap;
using std::pair;
using std::make_pair;
using std::less;
using std::string;
using std::type_info;
using std::exception;
using std::range_error;
using std::ios;
using std::ostream;
using std::istream;
using std::ofstream;
using std::ifstream;
using std::ostringstream;
using std::istringstream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::setprecision;
using std::setw;
using std::swap;
using std::min;
using std::max;
using std::mem_fun;
using std::sqrt;
//using std::pow;
using std::abs;
using std::atan2;
using std::isfinite;

/** Powers - standard or non-standard */
template <class ExponentT>
inline constexpr double pow(double x, ExponentT p) {
  return std::pow(x,double(p));
}

/** Square root of an integer. */
inline double sqrt(int x) {
  return std::sqrt(double(x));
}

/** factorial */
inline constexpr long double factorial(unsigned int n) {
  return (n < 2) ? 1.0 : n * factorial(n - 1);
}

/** Check if a given object is a part of a container. */
template <typename Container, typename Key>
inline bool member(const Container & c, const Key & k) {
  return c.find(k) != c.end();
}

/** Check if a given object is a part of a vector. */
template <typename T, typename Key>
inline bool member(const vector<T> & v, const Key & k) {
  for ( typename vector<T>::const_iterator i = v.begin(); i != v.end(); ++i )
    if ( *i == k ) return true;
  return false;
  // return find(v.begin(), v.end(), k) != v.end();
}

/** Return an insert iterator for a given container. */
template <typename Cont>
inline std::insert_iterator<Cont> inserter(Cont & c) {
  return std::insert_iterator<Cont>(c, c.end());
}


/** Return an insert iterator for a given vector. Overrides the
 *  general version. */
template <typename T, typename A>
inline std::back_insert_iterator< vector<T,A> > inserter(vector<T,A> & v) {
  return back_inserter(v);
}

/** Return an insert iterator for a given vector. Overrides the
 *  general version. */
template <typename T, typename A>
inline std::back_insert_iterator< deque<T,A> > inserter(deque<T,A> & v) {
  return back_inserter(v);
}

/** Stream manipulator setting an ostream to left-adjust its ouput. */
inline ostream& left(ostream& os) {
  os.setf(ios::left, ios::adjustfield);
  return os;
}

/** Stream manipulator setting an ostream to right-adjust its ouput. */
inline ostream& right(ostream& os) {
  os.setf(ios::right, ios::adjustfield);
  return os;
}

}

/** Macro for declaring a set. */
#define ThePEG_DECLARE_SET(VALTYPE,NAME)                               \
  /** A set of VALTYPE. */                                             \
  typedef set<VALTYPE, less<VALTYPE> > NAME

/** Macro for declaring a multiset. */
#define ThePEG_DECLARE_MULTISET(VALTYPE,NAME)                          \
  /** A multiset of VALTYPE. */                                        \
  typedef multiset<VALTYPE, less<VALTYPE> > NAME

/** Macro for declaring a map. */
#define ThePEG_DECLARE_MAP(KEYTYPE,VALTYPE,NAME)                       \
  /** A map of VALTYPE indexed by KEYTYPE. */                          \
  typedef map<KEYTYPE, VALTYPE, less<KEYTYPE> > NAME

#endif /* ThePEG_std_H */
