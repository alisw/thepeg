// -*- C++ -*-
//
// PolarizedBeamParticleData.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PolarizedBeamParticleData class.
//

#include "PolarizedBeamParticleData.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace ThePEG;

PolarizedBeamParticleData::PolarizedBeamParticleData(long newId, string newPDGName)
  : BeamParticleData(newId, newPDGName), theLongPolarization(0.0) {}

PDPtr PolarizedBeamParticleData::
Create(long newId, string newPDGName) {
  return new_ptr(PolarizedBeamParticleData(newId, newPDGName));
}

PDPair PolarizedBeamParticleData::
Create(long newId, string newPDGName, string newAntiPDGName) {
  PDPair pap;
  pap.first = new_ptr(PolarizedBeamParticleData(newId, newPDGName));
  pap.second = new_ptr(PolarizedBeamParticleData(-newId, newAntiPDGName));
  antiSetup(pap);
  return pap;
}

PDPtr PolarizedBeamParticleData::pdclone() const {
  return new_ptr(*this);
}

void PolarizedBeamParticleData::persistentOutput(PersistentOStream & os) const {
  os << theLongPolarization;
}

void PolarizedBeamParticleData::persistentInput(PersistentIStream & is, int) {
  is >> theLongPolarization;
}

ClassDescription<PolarizedBeamParticleData> 
PolarizedBeamParticleData::initPolarizedBeamParticleData;

void PolarizedBeamParticleData::Init() {

  static ClassDocumentation<PolarizedBeamParticleData> documentation
    ("There is no documentation for the ThePEG::PolarizedBeamParticleData class");

  static Parameter<PolarizedBeamParticleData,double> interfaceLongitudinalPolarization
    ("LongitudinalPolarization",
     "The longitudinal polarization",
     &PolarizedBeamParticleData::theLongPolarization, 0.0, -1.0, 1.0,
     false, false, Interface::limited);

  interfaceLongitudinalPolarization.rank(15);

}

RhoDMatrix PolarizedBeamParticleData::rhoMatrix() const {
  if(iSpin()!=PDT::Spin1Half) {
    throw Exception() << "Polarized Beams are currently only available for fermions\n"
		      << Exception::runerror;
  }
  RhoDMatrix output(PDT::Spin1Half);
  output(0,0) = 0.5*(1.-theLongPolarization);
  output(1,1) = 0.5*(1.+theLongPolarization);
  return output;
}
