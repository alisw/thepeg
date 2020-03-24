// -*- C++ -*-
//
// SamplerBase.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SamplerBase class.
//

#include "SamplerBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

SamplerBase::~SamplerBase() {}

void SamplerBase::persistentOutput(PersistentOStream & os) const {
  os << theEventHandler << theLastPoint;
  // Add all member variable which should be written persistently here.
}

void SamplerBase::persistentInput(PersistentIStream & is, int) {
  is >> theEventHandler >> theLastPoint;
  // Add all member variable which should be read persistently here.
}

AbstractClassDescription<SamplerBase> SamplerBase::initSamplerBase;
// Definition of the static class description member.

void SamplerBase::Init() {

  static ClassDocumentation<SamplerBase> documentation
    ("This is the base class for all phase space sampler classes to be"
     "used by the ThePEG::StandardEventHandler class to sample the phase"
     "space according to the cross sections for the proceses in the"
     "ThePEG::StandardEventHandler.");

}
