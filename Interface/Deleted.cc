// -*- C++ -*-
//
// Deleted.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DeletedBase class.
//

#include "InterfacedBase.h"
#include "Deleted.h"
#include "ThePEG/Utilities/Throw.h"

using namespace ThePEG;

string DeletedBase::exec(InterfacedBase &, string,
			 string) const {
  Throw<InterfaceException>()
    << "The interface '" << name() << "' has been removed. " << description()
    << Exception::runerror;
  return "";
}

string DeletedBase::type() const {
  return "Dd";
}

string DeletedBase::doxygenType() const {
  return "Deleted";
}

