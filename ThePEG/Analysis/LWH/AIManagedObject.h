// -*- C++ -*-
//
// AIManagedObject.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_AIManagedObject_H
#define LWH_AIManagedObject_H

#ifndef LWH_USING_AIDA

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

#else
#include "AIDA/IManagedObject.h"
#endif

#endif /* LWH_AIManagedObject_H */
