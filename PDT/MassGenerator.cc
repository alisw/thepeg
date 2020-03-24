// -*- C++ -*-
//
// MassGenerator.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MassGenerator class.
//

#include "MassGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

AbstractNoPIOClassDescription<MassGenerator> MassGenerator::initMassGenerator;

void MassGenerator::Init() {

  static ClassDocumentation<MassGenerator> documentation
    ("This is the base class for models giving specific masses to particle "
     "instances.");

}

