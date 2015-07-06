// -*- C++ -*-
//
// TypeInfo.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_TypeInfo_H
#define ThePEG_TypeInfo_H

#include "DescriptionList.h"

namespace ThePEG {

/** TypeInfo is a simple wrapper around the ClassDescription system in
 *  ThePEG. Defines a few static functions for returning the mane and
 *  version of the class given objects of that class. */
struct TypeInfo {

  /** Return the name of the class of the given object. */
  template <typename T>
  static string name(const T &)
  {
    const ClassDescriptionBase * cd = DescriptionList::find(typeid(T));
    if ( cd ) return cd->name();
    return "**** CLASS NOT REGISTERED ****";
  }

  /** Return the version number of the class of the given object. */
  template <typename T>
  static int version(const T &)
  {
    const ClassDescriptionBase * cd = DescriptionList::find(typeid(T));
    if ( cd ) return cd->version();
    return -1;
  }

};

}

#endif /* ThePEG_TypeInfo_H */
