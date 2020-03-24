// -*- C++ -*-
//
// DeltaMeasureCuts.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DeltaMeasureCuts class.
//

#include "DeltaMeasureCuts.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/PDT/MatcherBase.h"

using namespace ThePEG;

void DeltaMeasureCuts::describe() const {
  CurrentGenerator::log() 
    << fullName() << ":\n"
    << "MinDeltaR = " << theMinDeltaR << " \n"
    << "MinDeltaEta = " << theMinDeltaEta << " \n\n";
}

IBPtr DeltaMeasureCuts::clone() const {
  return new_ptr(*this);
}

IBPtr DeltaMeasureCuts::fullclone() const {
  return new_ptr(*this);
}

Energy DeltaMeasureCuts::minDeltaMeasureCuts(tcPDPtr, tcPDPtr) const {
  return ZERO;
}

Energy2 DeltaMeasureCuts::minSij(tcPDPtr, tcPDPtr) const {
  return ZERO;
}

Energy2 DeltaMeasureCuts::minTij(tcPDPtr, tcPDPtr) const {
  return ZERO;
}

Energy DeltaMeasureCuts::minKTClus(tcPDPtr, tcPDPtr) const {
  return ZERO;
}
double DeltaMeasureCuts::minDeltaR(tcPDPtr, tcPDPtr) const {
  return theMinDeltaR;
}

double DeltaMeasureCuts::minDurham(tcPDPtr, tcPDPtr) const {
  return 0.0;
}

bool DeltaMeasureCuts::passCuts(tcCutsPtr, tcPDPtr pitype, tcPDPtr pjtype,
		      LorentzMomentum pi, LorentzMomentum pj,
		      bool inci, bool incj) const {
  if ( theMatcher && !theMatcher->matches(*pitype) ) return true;
  if ( theMatcher && !theMatcher->matches(*pjtype) ) return true;
  if ( inci || incj ) return true;
  else {
    double deta2 = sqr(pi.eta() - pj.eta());
    if (abs(pi.eta() - pj.eta()) <= theMinDeltaEta) return false;
    double dphi = abs(pi.phi() - pj.phi());
    if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
    double dr = sqrt(deta2 + sqr(dphi));
    if ( dr <= theMinDeltaR ) return false;
  }
  return true;
}

void DeltaMeasureCuts::persistentOutput(PersistentOStream & os) const {
  os << theMinDeltaEta << theMinDeltaR << theMatcher;
}

void DeltaMeasureCuts::persistentInput(PersistentIStream & is, int) {
  is >> theMinDeltaEta >> theMinDeltaR >> theMatcher;
}

ClassDescription<DeltaMeasureCuts> DeltaMeasureCuts::initDeltaMeasureCuts;
// Definition of the static class description member.

void DeltaMeasureCuts::Init() {

  static ClassDocumentation<DeltaMeasureCuts> documentation
    ("This clas implements the cuts relevant for the "
     "\\f$\\Delta R\\f$-measure in the longitudinally invariant "
     "kt-algorithm. By default the cut is only applied to coloured "
     "particles, but optionally it may be applied to all particle types.");

   static Parameter<DeltaMeasureCuts,double> interfaceMinDeltaR
   ("MinDeltaR",
     "The minimum allowed legoplot distance ",
     &DeltaMeasureCuts::theMinDeltaR, 0.7, 0.0, 10.0,
     false, false, Interface::limited);

   static Parameter<DeltaMeasureCuts,double> interfaceMinDeltaEta
   ("MinDeltaEta",
     "The minimum allowed rapidity separation ",
     &DeltaMeasureCuts::theMinDeltaEta, 0.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Reference<DeltaMeasureCuts,MatcherBase> interfaceMatcher
    ("Matcher",
     "If non-null only particles matching this object will be affected "
     "by the cut.",
     &DeltaMeasureCuts::theMatcher, true, false, true, true, false);


}

