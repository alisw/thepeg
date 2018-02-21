// -*- C++ -*-
//
// BeamParticleData.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BeamParticleData class.
//

#include "BeamParticleData.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace ThePEG;

BeamParticleData::BeamParticleData(long newId, string newPDGName)
  : ParticleData(newId, newPDGName) {}

PDPtr BeamParticleData::
Create(long newId, string newPDGName) {
  return new_ptr(BeamParticleData(newId, newPDGName));
}

PDPair BeamParticleData::
Create(long newId, string newPDGName, string newAntiPDGName) {
  PDPair pap;
  pap.first = new_ptr(BeamParticleData(newId, newPDGName));
  pap.second = new_ptr(BeamParticleData(-newId, newAntiPDGName));
  antiSetup(pap);
  return pap;
}

PDPtr BeamParticleData::pdclone() const {
  return new_ptr(*this);
}

void BeamParticleData::persistentOutput(PersistentOStream & os) const {
  os << thePDF;
}

void BeamParticleData::persistentInput(PersistentIStream & is, int) {
  is >> thePDF;
}

ClassDescription<BeamParticleData> BeamParticleData::initBeamParticleData;

void BeamParticleData::setPDF(PDFPtr pdf) {
  if ( pdf && !pdf->canHandle(tcPDPtr(dynamic_cast<const ParticleData *>(this))) )
    throw BeamParticleWrongPDF(name(), pdf? pdf->name(): string("<NULL>"));
  thePDF = pdf;
}

void BeamParticleData::Init() {

  static ClassDocumentation<BeamParticleData> documentation
    ("There is no documentation for the ThePEG::BeamParticleData class");

  static Reference<BeamParticleData,PDFBase> interfacePDF
    ("PDF",
     "The parton densities for this beam particle.",
     &BeamParticleData::thePDF, false, false, true, true,
     &BeamParticleData::setPDF, 0, 0);

  interfacePDF.rank(15);

}

BeamParticleWrongPDF::BeamParticleWrongPDF(string p, string pdf) {
  theMessage << "The parton density object '" << pdf << "' cannot be used to "
	     << "handle densities of particle '" << p << "'. (Possibly due to "
	     << "the remnant handler assigned to the parton density.)";
  severity(warning);
}
