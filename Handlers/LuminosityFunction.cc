// -*- C++ -*-
//
// LuminosityFunction.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LuminosityFunction class.
//

#include "LuminosityFunction.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Vectors/LorentzRotation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "LuminosityFunction.tcc"
#endif

using namespace ThePEG;

IBPtr LuminosityFunction::clone() const {
  return new_ptr(*this);
}

IBPtr LuminosityFunction::fullclone() const {
  return new_ptr(*this);
}

LuminosityFunction::LuminosityFunction(Energy a, Energy b)
  : theBeamEMaxA(a), theBeamEMaxB(b) {}

void LuminosityFunction::select(tXCombPtr xcomb) {
  theLastXComb = xcomb;
}

bool LuminosityFunction::canHandle(const cPDPair &) const {
  return true;
}

Energy LuminosityFunction::maximumCMEnergy() const {
  return sqrt(4.0*beamEMaxA()*beamEMaxB());
}

LorentzRotation LuminosityFunction::getBoost() const {
  LorentzRotation r(0.0, 0.0, (beamEMaxA() - beamEMaxB())/
		              (beamEMaxA() + beamEMaxB()));
  return r;
}

double LuminosityFunction::Y() const {
  return 0.5*log(beamEMaxA()/beamEMaxB());
}

int LuminosityFunction::nDim(const cPDPair &) const {
  return 0;
}

double LuminosityFunction::
value(const cPDPair &, double l1, double l2) const {
  return l1 == 0.0 && l2 == 0.0? 1.0: 0.0;
}

pair<double,double>
LuminosityFunction::
generateLL(const double *, double &) const {
  return make_pair(0.0, 0.0);
}

void LuminosityFunction::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb << ounit(theBeamEMaxA, GeV) << ounit(theBeamEMaxB, GeV);
}

void LuminosityFunction::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb >> iunit(theBeamEMaxA, GeV) >> iunit(theBeamEMaxB, GeV);
}

ClassDescription<LuminosityFunction> LuminosityFunction::initLuminosityFunction;

void LuminosityFunction::Init() {

  static ClassDocumentation<LuminosityFunction> documentation
    ("This base class should be used by all classes describing the luminosity "
     "and energy distribution of colliding particle beams.");

  static Parameter<LuminosityFunction,Energy> interfaceBeamEMaxA
    ("BeamEMaxA",
     "The maximum energy of the beam entering along the positive z-axis. "
     "Note that derived classes may shift the beams away from the z-axis.",
     &LuminosityFunction::theBeamEMaxA, GeV, 45.6*GeV, ZERO, ZERO,
     true, false, Interface::lowerlim);

  static Parameter<LuminosityFunction,Energy> interfaceBeamEMaxB
    ("BeamEMaxB",
     "The maximum energy of the beam entering along the negative z-axis. "
     "Note that derived classes may shift the beams away from the z-axis.",
     &LuminosityFunction::theBeamEMaxB, GeV, 45.6*GeV, ZERO, ZERO,
     true, false, Interface::lowerlim);

  interfaceBeamEMaxA.rank(10);
  interfaceBeamEMaxB.rank(9);
  interfaceBeamEMaxA.setHasDefault(false);
  interfaceBeamEMaxB.setHasDefault(false);

}

