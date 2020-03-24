// -*- C++ -*-
//
// JetPairRegion.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
// Copyright (C) 2009-2019 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the JetPairRegion class.
//

#include "JetPairRegion.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/CurrentGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

JetPairRegion::JetPairRegion()
  : theMassMin(0.*GeV), theMassMax(Constants::MaxEnergy),
    theDeltaRMin(0.0), theDeltaRMax(Constants::MaxRapidity),
    theDeltaYMin(0.0), theDeltaYMax(Constants::MaxRapidity), 
    theOppositeHemispheres(false), theCutWeight(1.0) {}

JetPairRegion::~JetPairRegion() {}

IBPtr JetPairRegion::clone() const {
  return new_ptr(*this);
}

IBPtr JetPairRegion::fullclone() const {
  return new_ptr(*this);
}

void JetPairRegion::describe() const {

  CurrentGenerator::log()
    << "JetPairRegion '" << name() << "' matching "
    << " JetRegions '" << firstRegion()->name()
    << "' and '" << secondRegion()->name() << "' with\n";

  CurrentGenerator::log()
    << "m    = " << massMin()/GeV << " .. " << massMax()/GeV << " GeV\n"
    << "dR   = " << deltaRMin() << " .. " << deltaRMax() << "\n"
    << "dy   = " << deltaYMin() << " .. " << deltaYMax() << "\n";

}

bool JetPairRegion::matches(tcCutsPtr parent) {
  
  if ( !firstRegion()->didMatch() ||
       !secondRegion()->didMatch() ) {
    theCutWeight = 0.0;
    return false;
  }

  theCutWeight = 1.0;

  const LorentzMomentum& pi = firstRegion()->lastMomentum();
  const LorentzMomentum& pj = secondRegion()->lastMomentum();

  Energy m = (pi+pj).m();
  if ( !(firstRegion()->lessThanEnergy(massMin(),m,theCutWeight) &&
	 firstRegion()->lessThanEnergy(m,massMax(),theCutWeight)) )
    return false;

  double dy = abs(pi.rapidity() - pj.rapidity());
  if ( !(firstRegion()->lessThanRapidity(deltaYMin(),dy,theCutWeight) &&
	 firstRegion()->lessThanRapidity(dy,deltaYMax(),theCutWeight)) )
    return false;

  double dphi = abs(pi.phi() - pj.phi());
  if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
  double dR = sqrt(sqr(dy)+sqr(dphi));
  if ( !(firstRegion()->lessThanRapidity(deltaRMin(),dR,theCutWeight) &&
	 firstRegion()->lessThanRapidity(dR,deltaRMax(),theCutWeight)) )
    return false;

  double py = 
    (pi.rapidity() + parent->currentYHat()) * 
    (pj.rapidity() + parent->currentYHat());
  if ( theOppositeHemispheres && 
       !firstRegion()->lessThanRapidity(py,0.0,theCutWeight) )
    return false;

  return true;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void JetPairRegion::persistentOutput(PersistentOStream & os) const {
  os << theFirstRegion << theSecondRegion
     << ounit(theMassMin,GeV) << ounit(theMassMax,GeV)
     << theDeltaRMin << theDeltaRMax 
     << theDeltaYMin << theDeltaYMax 
     << theOppositeHemispheres << theCutWeight;
}

void JetPairRegion::persistentInput(PersistentIStream & is, int) {
  is >> theFirstRegion >> theSecondRegion
     >> iunit(theMassMin,GeV) >> iunit(theMassMax,GeV)
     >> theDeltaRMin >> theDeltaRMax 
     >> theDeltaYMin >> theDeltaYMax 
     >> theOppositeHemispheres >> theCutWeight;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<JetPairRegion,HandlerBase>
  describeThePEGJetPairRegion("ThePEG::JetPairRegion", "JetCuts.so");

void JetPairRegion::Init() {

  static ClassDocumentation<JetPairRegion> documentation
    ("JetPairRegion implements constraints on jets matching two jet regions.");

  
  static Reference<JetPairRegion,JetRegion> interfaceFirstRegion
    ("FirstRegion",
     "The first region to act on.",
     &JetPairRegion::theFirstRegion, false, false, true, false, false);

  static Reference<JetPairRegion,JetRegion> interfaceSecondRegion
    ("SecondRegion",
     "The second region to act on.",
     &JetPairRegion::theSecondRegion, false, false, true, false, false);

  static Parameter<JetPairRegion,Energy> interfaceMassMin
    ("MassMin",
     "The minimum jet-jet invariant mass.",
     &JetPairRegion::theMassMin, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<JetPairRegion,Energy> interfaceMassMax
    ("MassMax",
     "The maximum jet-jet invariant mass.",
     &JetPairRegion::theMassMax, GeV, Constants::MaxEnergy, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<JetPairRegion,double> interfaceDeltaRMin
    ("DeltaRMin",
     "The minimum jet-jet lego-plot separation.",
     &JetPairRegion::theDeltaRMin, 0.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<JetPairRegion,double> interfaceDeltaRMax
    ("DeltaRMax",
     "The maximum jet-jet lego-plot separation.",
     &JetPairRegion::theDeltaRMax, Constants::MaxRapidity, 0, 0,
     false, false, Interface::lowerlim);

  static Parameter<JetPairRegion,double> interfaceDeltaYMin
    ("DeltaYMin",
     "The minimum jet-jet rapidity separation.",
     &JetPairRegion::theDeltaYMin, 0.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<JetPairRegion,double> interfaceDeltaYMax
    ("DeltaYMax",
     "The maximum jet-jet rapidity separation.",
     &JetPairRegion::theDeltaYMax, Constants::MaxRapidity, 0, 0,
     false, false, Interface::lowerlim);

  static Switch<JetPairRegion,bool> interfaceOppositeHemispheres
    ("OppositeHemispheres",
     "Should the jets go into opposite detector hemispheres?",
     &JetPairRegion::theOppositeHemispheres, false, true, false);
  static SwitchOption interfaceOppositeHemispheresTrue
    (interfaceOppositeHemispheres,
     "True",
     "Require jets to be in different hemispheres",
     true);
  static SwitchOption interfaceOppositeHemispheresFalse
    (interfaceOppositeHemispheres,
     "False",
     "Do not require jets to be in different hemispheres",
     false);
  static SwitchOption interfaceOppositeHemispheresYes
    (interfaceOppositeHemispheres,
     "Yes",
     "Require jets to be in different hemispheres",
     true);
  static SwitchOption interfaceOppositeHemispheresNo
    (interfaceOppositeHemispheres,
     "No",
     "Do not require jets to be in different hemispheres",
     false);

}

