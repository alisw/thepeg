// -*- C++ -*-
//
// AlphaEMBase.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AlphaEMBase class.
//

#include "AlphaEMBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

AbstractNoPIOClassDescription<AlphaEMBase> AlphaEMBase::initAlphaEMBase;

void AlphaEMBase::Init() {

  static ClassDocumentation<AlphaEMBase> documentation
    ("An abstract base class used by the StandardModelBase class to "
     "implement the electro-magnetic coupling.");

}

