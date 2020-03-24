// -*- C++ -*-
//
// ThePEG.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

// This file contains the implementations of the declarations in
// ThePEG.h.

#include "ThePEG/Config/ThePEG.h"

using namespace ThePEG;

void Base::debug() const {
  debugme();
}

void Base::debugme() const {
  cerr << "(#ref: " << referenceCount() << ")";
}

