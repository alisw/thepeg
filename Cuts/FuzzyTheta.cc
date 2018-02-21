// -*- C++ -*-
//
// FuzzyTheta.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
// Copyright (C) 2009-2017 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FuzzyTheta class.
//

#include "FuzzyTheta.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

FuzzyTheta::FuzzyTheta() 
  : theEnergyWidth(1*GeV), theRapidityWidth(0.1),
    theAngularWidth(0.1) {}

FuzzyTheta::~FuzzyTheta() {}

IBPtr FuzzyTheta::clone() const {
  return new_ptr(*this);
}

IBPtr FuzzyTheta::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FuzzyTheta::persistentOutput(PersistentOStream & os) const {
  os << ounit(theEnergyWidth,GeV) << theRapidityWidth << theAngularWidth;
}

void FuzzyTheta::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theEnergyWidth,GeV) >> theRapidityWidth >> theAngularWidth;
}

ClassDescription<FuzzyTheta> FuzzyTheta::initFuzzyTheta;

void FuzzyTheta::Init() {

  static ClassDocumentation<FuzzyTheta> documentation
    ("FuzzyTheta implements fuzzy cut prescriptions.");

  static Parameter<FuzzyTheta,Energy> interfaceEnergyWidth
    ("EnergyWidth",
     "The width of smeared energy- and momentum-type cuts.",
     &FuzzyTheta::theEnergyWidth, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<FuzzyTheta,double> interfaceRapidityWidth
    ("RapidityWidth",
     "The width of smeared rapidity-type cuts.",
     &FuzzyTheta::theRapidityWidth, 0.1, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<FuzzyTheta,double> interfaceAngularWidth
    ("AngularWidth",
     "The width of smeared angular-type cuts.",
     &FuzzyTheta::theAngularWidth, 0.1, 0.0, Constants::pi,
     false, false, Interface::limited);

}

