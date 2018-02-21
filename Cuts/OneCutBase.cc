// -*- C++ -*-
//
// OneCutBase.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneCutBase class.
//

#include "OneCutBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace ThePEG;

OneCutBase::~OneCutBase() {}

void OneCutBase::describe() const {
  CurrentGenerator::log() << fullName() << " has no description.\n\n";
}

Energy OneCutBase::minMaxKT(tcPDPtr p) const {
  return minKT(p);
}

double OneCutBase::minMaxEta(tcPDPtr p) const {
  return minEta(p);
}

double OneCutBase::maxMinEta(tcPDPtr p) const {
  return maxEta(p);
}

bool OneCutBase::passCuts(tcCutsPtr parent,
			  tcPDPtr ptype, LorentzMomentum p) const {
  if ( p.perp() <= minKT(ptype) ) return false;
  double y = p.rapidity() + parent->Y() + parent->currentYHat();
  if ( p.mt()*sinh(y) <= p.perp()*sinh(minEta(ptype)) ) return false;
  if ( p.mt()*sinh(y) >= p.perp()*sinh(maxEta(ptype)) ) return false;
  return true;
}

bool OneCutBase::passCuts(tcCutsPtr parent, tcPPtr p) const {
  return passCuts(parent, p->dataPtr(), p->momentum());
}

Energy OneCutBase::minKT(tcPDPtr) const {
  return ZERO;
}

double OneCutBase::minEta(tcPDPtr) const {
  return -Constants::MaxRapidity;
}

double OneCutBase::maxEta(tcPDPtr) const {
  return Constants::MaxRapidity;
}

double OneCutBase::minRapidityMax(tcPDPtr) const {
  return -Constants::MaxRapidity;
}

double OneCutBase::maxRapidityMin(tcPDPtr) const {
  return Constants::MaxRapidity;
}

AbstractNoPIOClassDescription<OneCutBase> OneCutBase::initOneCutBase;
// Definition of the static class description member.

void OneCutBase::Init() {

  static ClassDocumentation<OneCutBase> documentation
    ("This class corresponds to a kinematical cut to be made on a single "
     "outgoing parton from a hard sub-process.");

}

