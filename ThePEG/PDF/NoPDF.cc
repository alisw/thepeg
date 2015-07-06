// -*- C++ -*-
//
// NoPDF.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NoPDF class.
//

#include "NoPDF.h"
#include "ThePEG/Utilities/Interval.h"
#include "ThePEG/PDF/RemnantHandler.h"
#include "ThePEG/PDT/ParticleData.h"

using namespace ThePEG;

IBPtr NoPDF::clone() const {
  return new_ptr(*this);
}

IBPtr NoPDF::fullclone() const {
  return new_ptr(*this);
}

bool NoPDF::canHandleParticle(tcPDPtr) const {
  return true;
}

bool NoPDF::canHandle(tcPDPtr particle) const {
  return canHandleParticle(particle) && remnantHandler() &&
    remnantHandler()->canHandle(particle, cPDVector());
}

bool NoPDF::hasPoleIn1(tcPDPtr particle, tcPDPtr parton) const {
  return particle == parton;
}

cPDVector NoPDF::partons(tcPDPtr p) const {
  return cPDVector(1, p);
}

double NoPDF::
xfl(tcPDPtr particle, tcPDPtr parton, Energy2, double l,
    Energy2) const {
  return ( l == 0 && particle == parton )? 1.0: 0.0;
}

NoPIOClassDescription<NoPDF> NoPDF::initNoPDF;

void NoPDF::Init() {}

