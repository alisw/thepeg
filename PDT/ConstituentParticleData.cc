// -*- C++ -*-
//
// ConstituentParticleData.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ConstituentParticleData class.
//

#include "ConstituentParticleData.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace ThePEG;

ConstituentParticleData::ConstituentParticleData(long newId, string newPDGName)
  : ParticleData(newId, newPDGName), theConstituentMass(ZERO),
    theDefaultConstituentMass(ZERO) {}


PDPtr ConstituentParticleData::
Create(long newId, string newPDGName) {
  return new_ptr(ConstituentParticleData(newId, newPDGName));
}

PDPair ConstituentParticleData::
Create(long newId, string newPDGName, string newAntiPDGName) {
  PDPair pap;
  pap.first = new_ptr(ConstituentParticleData(newId, newPDGName));
  pap.second = new_ptr(ConstituentParticleData(-newId, newAntiPDGName));
  antiSetup(pap);
  return pap;
}

void ConstituentParticleData::readSetup(istream & is) {
  ParticleData::readSetup(is);
  is >> iunit(theDefaultConstituentMass, GeV);
  theConstituentMass = theDefaultConstituentMass;
}

PDPtr ConstituentParticleData::pdclone() const {
  return new_ptr(*this);
}

void ConstituentParticleData::persistentOutput(PersistentOStream & os) const {
  os << ounit(theConstituentMass, GeV) << ounit(theDefaultConstituentMass, GeV);
}

void ConstituentParticleData::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theConstituentMass, GeV) >> iunit(theDefaultConstituentMass, GeV);
}

void ConstituentParticleData::setConstituentMass(Energy m) {
  theConstituentMass = m;
  ConstituentParticleData * apd =
    dynamic_cast<ConstituentParticleData*>(CC().operator->());
  if ( synchronized() && apd ) apd->theConstituentMass = m;
}

Energy ConstituentParticleData::defConstituentMass() const {
  return theDefaultConstituentMass;
}

ClassDescription<ConstituentParticleData>
ConstituentParticleData::initConstituentParticleData;

void ConstituentParticleData::Init() {

  static ClassDocumentation<ConstituentParticleData> documentation
    ("There is no documentation for the ThePEG::ConstituentParticleData class");

  static Parameter<ConstituentParticleData,Energy> interfaceMass
    ("ConstituentMass",
     "The constituent mass of the particle in GeV.",
     &ConstituentParticleData::theConstituentMass,
     GeV, ZERO, ZERO, Constants::MaxEnergy,
     false, false, Interface::lowerlim,
     &ConstituentParticleData::setConstituentMass, 0, 0, 0,
     &ConstituentParticleData::defConstituentMass);

  interfaceMass.rank(11.5);

}

