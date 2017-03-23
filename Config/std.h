// -*- C++ -*-
//
// std.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
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

/** @cond SHOWOSXWORKAROUNDS */
// Workarounds for OS X
#if defined __APPLE__ && defined __MACH__
extern "C" int isnan(double) throw();
extern "C" int isinf(double) throw();
#endif
/** @endcond */

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
using std::atan2;

/** Powers - standard or non-standard */
template <class ExponentT>
double pow(double x, ExponentT p) {
  return std::pow(x,double(p));
}

/** Square root of an integer. */
inline double sqrt(int x) {
  return std::sqrt(double(x));
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

#ifndef ThePEG_WRAP_STL_CONTAINERS

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

/** Macro for implementing a set. */
#define ThePEG_IMPLEMENT_SET(VALTYPE,NAME)

/** Macro for implementing a multiset. */
#define ThePEG_IMPLEMENT_MULTISET(VALTYPE,NAME)

/** Macro for implementing a map. */
#define ThePEG_IMPLEMENT_MAP(KEYTYPE,VALTYPE,NAME)

#else

/** Macro for declaring a set. */
#define ThePEG_DECLARE_SET(VALTYPE,NAME)                                \
class NAME : public set<VALTYPE, less<VALTYPE> > {  \
public:                                                                 \
  typedef set<VALTYPE, less<VALTYPE> > SETTYPE;     \
  NAME();                                                               \
  explicit NAME(const key_compare & c,                                  \
		const allocator_type & a = allocator_type());           \
  template <typename InputIterator>                                     \
  NAME(InputIterator first, InputIterator last)                         \
    : SETTYPE(first, last) {}                                           \
  template <typename InputIterator>                                     \
  NAME(InputIterator first, InputIterator last, const key_compare & c,  \
       const allocator_type & a = allocator_type())                     \
    : SETTYPE(first, last, c, a) {}                                     \
  NAME(const SETTYPE & s);                                              \
  NAME(const NAME & s);                                                 \
  ~NAME();                                                              \
  NAME & operator=(const NAME &);                                       \
  NAME & operator=(const SETTYPE &);                                    \
  pair<iterator,bool> insert(const value_type & x);                     \
  iterator insert(iterator position, const value_type & x);             \
  template <typename InputIterator>                                     \
  void insert(InputIterator first, InputIterator last) {                \
    SETTYPE::insert(first, last);                                       \
  }                                                                     \
  void erase(iterator position);                                        \
  size_type erase(const key_type & x);                                  \
  void erase(iterator first, iterator last);                            \
  void clear();                                                         \
  iterator find(const key_type & x) const;                              \
  size_type count(const key_type & x) const;                            \
  iterator lower_bound(const key_type & x) const;                       \
  iterator upper_bound(const key_type & x) const;                       \
  pair<iterator,iterator> equal_range(const key_type & x) const;        \
}


/** Macro for declaring a multiset. */
#define ThePEG_DECLARE_MULTISET(VALTYPE,NAME)                           \
class NAME :                                                            \
 public multiset<VALTYPE, less<VALTYPE> > {         \
public:                                                                 \
  typedef multiset<VALTYPE, less<VALTYPE> > SETTYPE;\
  NAME();                                                               \
  explicit NAME(const key_compare & c,                                  \
		const allocator_type & a = allocator_type());           \
  template <typename InputIterator>                                     \
  NAME(InputIterator first, InputIterator last)                         \
    : SETTYPE(first, last) {}                                           \
  template <typename InputIterator>                                     \
  NAME(InputIterator first, InputIterator last, const key_compare & c,  \
       const allocator_type & a = allocator_type())                     \
    : SETTYPE(first, last, c, a) {}                                     \
  NAME(const SETTYPE & s);                                              \
  NAME(const NAME & s);                                                 \
  ~NAME();                                                              \
  NAME & operator=(const NAME &);                                       \
  NAME & operator=(const SETTYPE &);                                    \
  iterator insert(const value_type & x);                                \
  iterator insert(iterator position, const value_type & x);             \
  template <typename InputIterator>                                     \
  void insert(InputIterator first, InputIterator last) {                \
    SETTYPE::insert(first, last);                                       \
  }                                                                     \
  void erase(iterator position);                                        \
  size_type erase(const key_type & x);                                  \
  void erase(iterator first, iterator last);                            \
  void clear();                                                         \
  iterator find(const key_type & x) const;                              \
  size_type count(const key_type & x) const;                            \
  iterator lower_bound(const key_type & x) const;                       \
  iterator upper_bound(const key_type & x) const;                       \
  pair<iterator,iterator> equal_range(const key_type & x) const;        \
}


/** Macro for declaring a map. */
#define ThePEG_DECLARE_MAP(KEYTYPE,VALTYPE,NAME)                        \
class NAME :                                                            \
  public map<KEYTYPE, VALTYPE, less<KEYTYPE> > {                        \
public:                                                                 \
  typedef map<KEYTYPE, VALTYPE, less<KEYTYPE> > MAPTYPE;		\
  NAME();                                                               \
  explicit NAME(const key_compare & c,                                  \
		const allocator_type & a = allocator_type());           \
  template <typename InputIterator>                                     \
  NAME(InputIterator first, InputIterator last)                         \
    : MAPTYPE(first, last) {}                                           \
  template <typename InputIterator>                                     \
  NAME(InputIterator first, InputIterator last, const key_compare & c,  \
       const allocator_type & a = allocator_type())                     \
    : MAPTYPE(first, last, c, a) {}                                     \
  NAME(const NAME & s);                                                 \
  NAME(const MAPTYPE & s);                                              \
  ~NAME();                                                              \
  NAME & operator=(const NAME &);                                       \
  NAME & operator=(const MAPTYPE &);                                    \
  data_type & operator[](const key_type & k);                           \
  pair<iterator,bool> insert(const value_type & x);                     \
  iterator insert(iterator position, const value_type & x);             \
  template <typename InputIterator>                                     \
  void insert(InputIterator first, InputIterator last) {                \
    MAPTYPE::insert(first, last);                                       \
  }                                                                     \
  void erase(iterator position);                                        \
  size_type erase(const key_type & x);                                  \
  void erase(iterator first, iterator last);                            \
  void clear();                                                         \
  iterator find(const key_type & x);                                    \
  const_iterator find(const key_type & x) const;                        \
  size_type count(const key_type & x) const;                            \
  iterator lower_bound(const key_type & x);                             \
  const_iterator lower_bound(const key_type & x) const;                 \
  iterator upper_bound(const key_type & x);                             \
  const_iterator upper_bound(const key_type & x) const;                 \
  pair<iterator,iterator> equal_range(const key_type & x);              \
  pair<const_iterator,const_iterator>                                   \
  equal_range(const key_type & x) const;                                \
}


/** Macro for implementing a set. */
#define ThePEG_IMPLEMENT_SET(VALTYPE,NAME)                             \
NAME::NAME() {}                                                         \
NAME::NAME(const key_compare & c, const allocator_type & a)             \
  :SETTYPE(c, a) {}                                                     \
NAME::NAME(const NAME & x) : SETTYPE(x) {}                              \
NAME::NAME(const SETTYPE & x) : SETTYPE(x) {}                           \
NAME::~NAME() {}                                                        \
NAME & NAME::operator=(const NAME & x) {                                \
  SETTYPE::operator=(x);                                                \
  return *this;                                                         \
}                                                                       \
NAME & NAME::operator=(const SETTYPE & x) {                             \
  SETTYPE::operator=(x);                                                \
 return *this;                                                          \
}                                                                       \
pair<NAME::iterator,bool> NAME::insert(const value_type & x) {          \
  return SETTYPE::insert(x);                                            \
}                                                                       \
NAME::iterator NAME::insert(iterator position, const value_type & x) {  \
  return SETTYPE::insert(position, x);                                  \
}                                                                       \
void NAME::erase(iterator position) {                                   \
  SETTYPE::erase(position);                                             \
}                                                                       \
NAME::size_type NAME::erase(const key_type & x) {                       \
  return SETTYPE::erase(x);                                             \
}                                                                       \
void NAME::erase(iterator first, iterator last) {                       \
  SETTYPE::erase(first, last);                                          \
}                                                                       \
void NAME::clear() {                                                    \
  SETTYPE::clear();                                                     \
}                                                                       \
NAME::iterator NAME::find(const key_type & x) const {                   \
  return SETTYPE::find(x);                                              \
}                                                                       \
NAME::size_type NAME::count(const key_type & x) const {                 \
  return SETTYPE::count(x);                                             \
}                                                                       \
NAME::iterator NAME::lower_bound(const key_type & x) const {            \
  return SETTYPE::lower_bound(x);                                       \
}                                                                       \
NAME::iterator NAME::upper_bound(const key_type & x) const {            \
  return SETTYPE::upper_bound(x);                                       \
}                                                                       \
pair<NAME::iterator,NAME::iterator>                                     \
NAME::equal_range(const key_type & x) const {                           \
  return SETTYPE::equal_range(x);                                       \
}                                                                       \


/** Macro for implementing a multiset. */
#define ThePEG_IMPLEMENT_MULTISET(VALTYPE,NAME)                        \
NAME::NAME() {}                                                         \
NAME::NAME(const key_compare & c, const allocator_type & a)             \
  :SETTYPE(c, a) {}                                                     \
NAME::NAME(const NAME & x) : SETTYPE(x) {}                              \
NAME::NAME(const SETTYPE & x) : SETTYPE(x) {}                           \
NAME::~NAME() {}                                                        \
NAME & NAME::operator=(const NAME & x) {                                \
  SETTYPE::operator=(x);                                                \
  return *this;                                                         \
}                                                                       \
NAME & NAME::operator=(const SETTYPE & x) {                             \
  SETTYPE::operator=(x);                                                \
 return *this;                                                          \
}                                                                       \
NAME::iterator NAME::insert(const value_type & x) {                     \
  return SETTYPE::insert(x);                                            \
}                                                                       \
NAME::iterator NAME::insert(iterator position, const value_type & x) {  \
  return SETTYPE::insert(position, x);                                  \
}                                                                       \
void NAME::erase(iterator position) {                                   \
  SETTYPE::erase(position);                                             \
}                                                                       \
NAME::size_type NAME::erase(const key_type & x) {                       \
  return SETTYPE::erase(x);                                             \
}                                                                       \
void NAME::erase(iterator first, iterator last) {                       \
  SETTYPE::erase(first, last);                                          \
}                                                                       \
void NAME::clear() {                                                    \
  SETTYPE::clear();                                                     \
}                                                                       \
NAME::iterator NAME::find(const key_type & x) const {                   \
  return SETTYPE::find(x);                                              \
}                                                                       \
NAME::size_type NAME::count(const key_type & x) const {                 \
  return SETTYPE::count(x);                                             \
}                                                                       \
NAME::iterator NAME::lower_bound(const key_type & x) const {            \
  return SETTYPE::lower_bound(x);                                       \
}                                                                       \
NAME::iterator NAME::upper_bound(const key_type & x) const {            \
  return SETTYPE::upper_bound(x);                                       \
}                                                                       \
pair<NAME::iterator,NAME::iterator>                                     \
NAME::equal_range(const key_type & x) const {                           \
  return SETTYPE::equal_range(x);                                       \
}                                                                       \


/** Macro for implementing a map. */
#define ThePEG_IMPLEMENT_MAP(KEYTYPE,VALTYPE,NAME)                     \
NAME::NAME() {}                                                         \
NAME::NAME(const key_compare & c, const allocator_type & a)             \
  :MAPTYPE(c, a) {}                                                     \
NAME::NAME(const NAME & x) : MAPTYPE(x) {}                              \
NAME::NAME(const MAPTYPE & x) : MAPTYPE(x) {}                            \
NAME::~NAME() {}                                                        \
NAME & NAME::operator=(const NAME & x) {                                \
  MAPTYPE::operator=(x);                                                \
  return *this;                                                         \
}                                                                       \
NAME & NAME::operator=(const MAPTYPE & x) {                             \
  MAPTYPE::operator=(x);                                                \
 return *this;                                                          \
}                                                                       \
pair<NAME::iterator,bool> NAME::insert(const value_type & x) {          \
  return MAPTYPE::insert(x);                                            \
}                                                                       \
NAME::iterator NAME::insert(iterator position, const value_type & x) {  \
  return MAPTYPE::insert(position, x);                                  \
}                                                                       \
void NAME::erase(iterator position) {                                   \
  MAPTYPE::erase(position);                                             \
}                                                                       \
NAME::size_type NAME::erase(const key_type & x) {                       \
  return MAPTYPE::erase(x);                                             \
}                                                                       \
void NAME::erase(iterator first, iterator last) {                       \
  MAPTYPE::erase(first, last);                                          \
}                                                                       \
void NAME::clear() {                                                    \
  MAPTYPE::clear();                                                     \
}                                                                       \
NAME::iterator NAME::find(const key_type & x) {                         \
  return MAPTYPE::find(x);                                              \
}                                                                       \
NAME::const_iterator NAME::find(const key_type & x) const {             \
  return MAPTYPE::find(x);                                              \
}                                                                       \
NAME::size_type NAME::count(const key_type & x) const {                 \
  return MAPTYPE::count(x);                                             \
}                                                                       \
NAME::iterator NAME::lower_bound(const key_type & x) {                  \
  return MAPTYPE::lower_bound(x);                                       \
}                                                                       \
NAME::const_iterator NAME::lower_bound(const key_type & x) const {      \
  return MAPTYPE::lower_bound(x);                                       \
}                                                                       \
NAME::iterator NAME::upper_bound(const key_type & x) {                  \
  return MAPTYPE::upper_bound(x);                                       \
}                                                                       \
NAME::const_iterator NAME::upper_bound(const key_type & x) const {      \
  return MAPTYPE::upper_bound(x);                                       \
}                                                                       \
pair<NAME::iterator,NAME::iterator>                                     \
NAME::equal_range(const key_type & x) {                                 \
  return MAPTYPE::equal_range(x);                                       \
}                                                                       \
pair<NAME::const_iterator,NAME::const_iterator>                         \
NAME::equal_range(const key_type & x) const {                           \
  return MAPTYPE::equal_range(x);                                       \
}                                                                       \
NAME::data_type & NAME::operator[](const key_type & k) {                \
  return MAPTYPE::operator[](k);                                        \
}                                                                       \

#endif

// #include "std.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "std.tcc"
#endif

#endif /* ThePEG_std_H */
