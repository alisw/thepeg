// -*- C++ -*-
//
// ManagedObject.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef LWH_ManagedObject_H
#define LWH_ManagedObject_H
//
// This is the declaration of the ManagedObject class.
//

#include "AIManagedObject.h"
#include <iostream>

namespace LWH {

using namespace AIDA;

/**
 * The creator of trees.
 */
class ManagedObject: public IManagedObject {

public:

  /// Destructor.
  virtual ~ManagedObject() {}

  /**
   * Write out the object to the given stream in XML format.
   */
  virtual bool writeXML(std::ostream & os,
			std::string path, std::string name) = 0;

  /**
   * Write out the object to the given stream in simple table format.
   */
  virtual bool writeFLAT(std::ostream & os,
			 std::string path, std::string name) = 0;

};

}

#endif /* LWH_ManagedObject_H */
