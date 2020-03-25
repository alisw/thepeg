// -*- C++ -*-
//
// AIManagedObject.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_AIManagedObject_H
#define LWH_AIManagedObject_H



#include <string>

/** @cond DONT_DOCUMENT_STRIPPED_DOWN_AIDA_INTERFACES */

namespace AIDA {

class ITree;

class IManagedObject {

public:

  virtual ~IManagedObject() {}

  virtual std::string name() const = 0;
  virtual void * cast(const std::string & className) const = 0;

};

}

/** @endcond */





#endif /* LWH_AIManagedObject_H */
