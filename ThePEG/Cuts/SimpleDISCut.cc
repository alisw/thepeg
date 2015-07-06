// -*- C++ -*-
//
// SimpleDISCut.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleDISCut class.
//

#include "SimpleDISCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace ThePEG;

void SimpleDISCut::describe() const {
  CurrentGenerator::log() 
    << fullName() << ":\n"
    << "Q2 = " << theMinQ2/GeV2 << " .. " << theMaxQ2/GeV2 << " GeV2\n"
    << "W2 = " << theMinW2/GeV2 << " .. " << theMaxW2/GeV2 << " GeV2\n\n";
}

IBPtr SimpleDISCut::clone() const {
  return new_ptr(*this);
}

IBPtr SimpleDISCut::fullclone() const {
  return new_ptr(*this);
}

bool SimpleDISCut::check(long idi, long ido) const {
  if ( abs(idi) <= 10 || abs(idi) > 16 ) return false;
  if ( idi*ido <= 0 ) return false;
  if ( chargedCurrent ) {
    if ( abs(idi)%2 ) return abs(idi) == abs(ido) - 1;
    else return abs(idi) == abs(ido) + 1;
  } else {
    return ( idi == ido );
  }
}

Energy2 SimpleDISCut::minTij(tcPDPtr pi, tcPDPtr po) const {
  if ( !check(pi->id(), po->id()) ) return ZERO;
  return theMinQ2;
}
  
bool SimpleDISCut::passCuts(tcCutsPtr cuts, tcPDPtr pitype, tcPDPtr pjtype,
			    LorentzMomentum pi, LorentzMomentum pj,
			    bool inci, bool incj) const {
  if ( inci ) {
    if ( incj ) return true;
    if ( !check(pitype->id(), pjtype->id()) ) return true;
    Energy2 Q2 = -(pi - pj).m2();
    double x = min(1.0, sqrt(cuts->currentSHat()/cuts->SMax())*
		   exp(-cuts->currentYHat()));
    Energy2 W2 = (1.0 - x)*Q2/x;
    return Q2 > theMinQ2 && Q2 < theMaxQ2 && W2 > theMinW2 && W2 < theMaxW2;
  }
  else if ( incj ) {
    if ( !check(pjtype->id(), pitype->id()) ) return true;
    Energy2 Q2 = -(pj - pi).m2();
    double x =
      min(1.0, sqrt(cuts->currentSHat()/cuts->SMax())*
	  exp(cuts->currentYHat()));
    Energy2 W2 = (1.0 - x)*Q2/x;
    return Q2 > theMinQ2 && Q2 < theMaxQ2 && W2 > theMinW2 && W2 < theMaxW2;
  }
  return true;
}

Energy2 SimpleDISCut::minSij(tcPDPtr, tcPDPtr) const {
  return ZERO;
}

double SimpleDISCut::minDeltaR(tcPDPtr, tcPDPtr) const {
  return 0.0;
}

Energy SimpleDISCut::minKTClus(tcPDPtr, tcPDPtr) const {
  return ZERO;
}

double SimpleDISCut::minDurham(tcPDPtr, tcPDPtr) const {
  return 0.0;
}

void SimpleDISCut::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMinQ2, GeV2) << ounit(theMaxQ2, GeV2)
     << ounit(theMinW2, GeV2) << ounit(theMaxW2, GeV2) << chargedCurrent;
}

void SimpleDISCut::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMinQ2, GeV2) >> iunit(theMaxQ2, GeV2)
     >> iunit(theMinW2, GeV2) >> iunit(theMaxW2, GeV2) >> chargedCurrent;
}

ClassDescription<SimpleDISCut> SimpleDISCut::initSimpleDISCut;
// Definition of the static class description member.

Energy2 SimpleDISCut::maxMinQ2() const {
  return theMaxQ2;
}

Energy2 SimpleDISCut::minMaxQ2() const {
  return theMinQ2;
}

Energy2 SimpleDISCut::maxMinW2() const {
  return theMaxW2;
}

Energy2 SimpleDISCut::minMaxW2() const {
  return theMinW2;
}

void SimpleDISCut::Init() {

  static ClassDocumentation<SimpleDISCut> documentation
    ("SimpleDISCut inherits from TwoCutBase and omplements a simple "
     "\\f$Q^2\\f$ cut on the a scattered lepton, either neutral or charged "
     "current.");


  static Parameter<SimpleDISCut,Energy2> interfaceMinQ2
    ("MinQ2",
     "The minimum \\f$Q^2\\f$.",
     &SimpleDISCut::theMinQ2, GeV2, 1.0*GeV2, ZERO, Constants::MaxEnergy2,
     true, false, Interface::limited,
     0, 0, 0, &SimpleDISCut::maxMinQ2, 0);

  static Parameter<SimpleDISCut,Energy2> interfaceMaxQ2
    ("MaxQ2",
     "The maximum \\f$Q^2\\f$. Note that this is only applied as a post-cut "
     "and will not affect the initial phase space cuts in the generation.",
     &SimpleDISCut::theMaxQ2, GeV2, 100.0*GeV2, ZERO, ZERO,
     true, false, Interface::lowerlim,
     0, 0, &SimpleDISCut::minMaxQ2, 0, 0);

  static Parameter<SimpleDISCut,Energy2> interfaceMinW2
    ("MinW2",
     "The minimum \\f$W^2\\f$. Note that this is only applied as a post-cut "
     "and will not affect the initial phase space cuts in the generation.",
     &SimpleDISCut::theMinW2, GeV2, 100.0*GeV2, ZERO, Constants::MaxEnergy2,
     true, false, Interface::limited,
     0, 0, 0, &SimpleDISCut::maxMinW2, 0);

  static Parameter<SimpleDISCut,Energy2> interfaceMaxW2
    ("MaxW2",
     "The maximum \\f$W^2\\f$. Note that this is only applied as a post-cut "
     "and will not affect the initial phase space cuts in the generation.",
     &SimpleDISCut::theMaxW2, GeV2, 1000000.0*GeV2, ZERO, ZERO,
     true, false, Interface::lowerlim,
     0, 0, &SimpleDISCut::minMaxW2, 0, 0);

  static Switch<SimpleDISCut,bool> interfaceCurrent
    ("Current",
     "Determines whether this cut should be applied to charged or neutral "
     "current events.",
     &SimpleDISCut::chargedCurrent, false, true, false);
  static SwitchOption interfaceCurrentCharged
    (interfaceCurrent,
     "Charged",
     "The cut is only applied to charged current events.",
     true);
  static SwitchOption interfaceCurrentNeutral
    (interfaceCurrent,
     "Neutral",
     "The cut is only applied to neutral current events.",
     false);

  interfaceMinQ2.rank(10);
  interfaceMaxQ2.rank(9);
  interfaceCurrent.rank(8);
  interfaceMinQ2.setHasDefault(false);
  interfaceMaxQ2.setHasDefault(false);
  interfaceMinW2.setHasDefault(false);
  interfaceMaxW2.setHasDefault(false);
  interfaceCurrent.setHasDefault(false);
}

