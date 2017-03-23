// -*- C++ -*-
//
// EnumIO.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_EnumIO_H
#define ThePEG_EnumIO_H
// This is the declaration of the IEnum and OEnum classes and
// associated templated functions.

// #include "ThePEG/Config/ThePEG.h"
// #include "EnumIO.fh"
// #include "EnumIO.xh"

namespace ThePEG {

template <typename T>
/**
 * The <code>OEnum</code> helper class is used to facilitate output of enums
 * to persistent streams. An enum can hence be written like this:<BR>
 * <code>os >> oenum(x);</code><BR>
 *
 * @see PersistentOStream
 * @see PersistentIStream
 * 
 */
struct OEnum {

  /** Constructor. */
  OEnum(const T & t): theT(t) {}

  /** Copy constructor. */
  OEnum(const OEnum & oe): theT(oe.theT) {}

  /** The variable to be written */
  const T & theT;

};

/**
 * The <code>IEnum</code> helper class is used to facilitate input of enums
 * from persistent streams. An enum can hence be read like this:<BR>
 * <code>is >> ienum(x);</code>
 *
 * @see PersistentOStream
 * @see PersistentIStream
 * 
 */
template <typename T>
struct IEnum {

  /** Constructor. */
  IEnum(T & t): theT(t) {}

  /** Copy constructor. */
  IEnum(const IEnum & ie): theT(ie.theT) {}

  /** The variable to be read */
  T & theT;

};

/** Helper function to create an OEnum object for a given variable. */
template <typename T>
inline OEnum<T> oenum(const T & t) {
  return OEnum<T>(t);
}

/** Helper function to create an IEnum object for a given variable. */
template <typename T>
inline IEnum<T> ienum(T & t) {
  return IEnum<T>(t);
}

/** Overloading of operator<< for OEnum. */
template <typename OStream, typename T>
OStream & operator<<(OStream & os, const OEnum<T> & e) {
  os << long(e.theT);
  return os;
}

/** Overloading of operator<< for IEnum. */
template <typename IStream, typename T>
IStream & operator>>(IStream & is, const IEnum<T> & e) {
  long l = 0;
  is >> l;
  e.theT = T(l);
  return is;
}

}

#endif /* ThePEG_EnumIO_H */
