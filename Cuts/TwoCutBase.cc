// -*- C++ -*-
//
// TwoCutBase.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoCutBase class.
//

#include "TwoCutBase.h"
#include "Cuts.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace ThePEG;

TwoCutBase::~TwoCutBase() {}

void TwoCutBase::describe() const {
  CurrentGenerator::log() << fullName() << " has no description.\n\n";
}

bool TwoCutBase::passCuts(tcCutsPtr parent, tcPDPtr pitype, tcPDPtr pjtype,
			  LorentzMomentum pi, LorentzMomentum pj,
			  bool inci, bool incj) const {
  if ( inci && incj ) return true;
  else if ( inci ) {
    if ( -(pj - pi).m2() <= minTij(pitype, pjtype) ) return false;
    if ( pj.perp() <= minKTClus(tcPDPtr(), pjtype) ) return false;
  }
  else if ( incj ) {
    if ( -(pi - pj).m2() <= minTij(pjtype, pitype) ) return false;
    if ( pi.perp() <= minKTClus(tcPDPtr(), pitype) ) return false;
  }
  else {
    if ( (pi + pj).m2() <= minSij(pitype, pjtype) ) return false;
    double deta2 = sqr(pi.eta() - pj.eta());
    double dphi = abs(pi.phi() - pj.phi());
    if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
    double dr = sqrt(deta2 + sqr(dphi));
    if ( dr < minDeltaR(pitype, pjtype) ) return false;
    if ( min(pi.perp(), pj.perp())*dr <= minKTClus(pitype, pjtype) )
      return false;
    if ( 2.0*sqr(min(pi.e(), pj.e()))*(1.0 - cos(pi.angle(pj))) <
	 parent->currentSHat()*minDurham(pitype, pjtype) ) return false;
  }
  return true;
}

bool TwoCutBase::passCuts(tcCutsPtr parent, tcPPtr pi, tcPPtr pj,
			  bool inci, bool incj) const {
  return passCuts(parent, pi->dataPtr(), pj->dataPtr(),
		  pi->momentum(), pj->momentum(), inci, incj);
}

AbstractNoPIOClassDescription<TwoCutBase> TwoCutBase::initTwoCutBase;
// Definition of the static class description member.

void TwoCutBase::Init() {

  static ClassDocumentation<TwoCutBase> documentation
    ("This class corresponds to a kinematical cut to be made on a pair of "
     "particles in a hard sub-process.");

}

