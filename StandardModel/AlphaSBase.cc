// -*- C++ -*-
//
// AlphaSBase.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AlphaSBase class.
//

#include "AlphaSBase.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

using namespace ThePEG;

void AlphaSBase::doinit() {
  theFlavourThresholds = flavourThresholds();
  theLambdaQCDs = LambdaQCDs();
  RunningCoupling::doinit();
}

void AlphaSBase::persistentOutput(PersistentOStream & os) const {
  os << ounit(theQuarkMasses,GeV) 
     << ounit(theFlavourThresholds, GeV2) << ounit(theLambdaQCDs, GeV);
}

void AlphaSBase::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theQuarkMasses,GeV) 
     >> iunit(theFlavourThresholds, GeV2) >> iunit(theLambdaQCDs, GeV);
}

AbstractClassDescription<AlphaSBase> AlphaSBase::initAlphaSBase;

void AlphaSBase::Init() {

  static ClassDocumentation<AlphaSBase> documentation
    ("An abstract base class used by the StandardModelBase to implement the "
     "QCD coupling.");


  static ParVector<AlphaSBase,Energy> interfaceQuarkMasses
    ("QuarkMasses",
     "The quark masses to be used instead of the masses set in the particle data.",
     &AlphaSBase::theQuarkMasses, GeV, -1, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

}

