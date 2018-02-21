// -*- C++ -*-
//
// ZGenerator.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//$Id$
// ----------------------------------------------------
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZGenerator class.
//

#include "ZGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

AbstractClassDescription<ZGenerator> ZGenerator::initZGenerator;

void ZGenerator::Init(){

  static ClassDocumentation<ZGenerator> documentation
    ("The base class for all classes implementing models to generate the "
     "momentum fraction, \\f$z\\f$, taken by hadrons produced in a "
     "hadronization scenario.");

}

