// -*- C++ -*-
//
// SimpleKTCut.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleKTCut class.
//

#include "SimpleKTCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace ThePEG;

SimpleKTCut::~SimpleKTCut() {}

void SimpleKTCut::describe() const {
  CurrentGenerator::log() 
    << fullName() << ":\n"
    << "KT  = " << theMinKT/GeV << " .. " << theMaxKT/GeV << " GeV\n"
    << "Eta = " << theMinEta << " .. " << theMaxEta << "\n\n";
}

IBPtr SimpleKTCut::clone() const {
  return new_ptr(*this);
}

IBPtr SimpleKTCut::fullclone() const {
  return new_ptr(*this);
}

Energy SimpleKTCut::minKT(tcPDPtr p) const {
  if ( theMatcher && !theMatcher->matches(*p) ) return ZERO;
  return theMinKT;
}

double SimpleKTCut::minEta(tcPDPtr p) const {
  if ( theMatcher && !theMatcher->matches(*p) )
    return -Constants::MaxRapidity;
  return theMinEta;
}

double SimpleKTCut::maxEta(tcPDPtr p) const {
  if ( theMatcher && !theMatcher->matches(*p) )
    return Constants::MaxRapidity;
  return theMaxEta;
}

bool SimpleKTCut::passCuts(tcCutsPtr parent,
			   tcPDPtr ptype, LorentzMomentum p) const {
  if ( theMatcher && !theMatcher->matches(*ptype) ) return true;
  if ( p.perp() < theMinKT ) return false;
  if ( p.perp() > theMaxKT ) return false;
  double y = abs(p.t()) <= abs(p.z()) ? (p .z() > ZERO ? 1e10 : -1e10) : p.rapidity();
  y += parent->Y() + parent->currentYHat();
  if ( p.mt()*sinh(y) <= p.perp()*sinh(theMinEta) ) return false;
  if ( p.mt()*sinh(y) >= p.perp()*sinh(theMaxEta) ) return false;
  return true;
}

void SimpleKTCut::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMinKT, GeV) << ounit(theMaxKT, GeV)
     << theMinEta << theMaxEta << theMatcher;
}

void SimpleKTCut::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMinKT, GeV) >> iunit(theMaxKT, GeV)
     >> theMinEta >> theMaxEta >> theMatcher;
}

ClassDescription<SimpleKTCut> SimpleKTCut::initSimpleKTCut;
// Definition of the static class description member.

Energy SimpleKTCut::maxKTMin() const {
  return theMaxKT;
}

Energy SimpleKTCut::minKTMax() const {
  return theMinKT;
}

double SimpleKTCut::maxEtaMin() const {
  return theMaxEta;
}

double SimpleKTCut::minEtaMax() const {
  return theMinEta;
}

void SimpleKTCut::Init() {

  typedef double (ThePEG::SimpleKTCut::*IGFN)() const;
  typedef void (ThePEG::SimpleKTCut::*ISFN)(double);
  typedef Energy (ThePEG::SimpleKTCut::*IGFNK)() const;
  typedef void (ThePEG::SimpleKTCut::*ISFNK)(Energy);

  static ClassDocumentation<SimpleKTCut> documentation
    ("This is a very simple concrete sub-class of OneCutbase simply "
     "requiring a minimum transverse momentum of any outgoing particle. "
     "It is also possible to require a minimum and maximum pseudorapidity. "
     "Optionally the restrictions only apply to particles matching a "
     "specific matcher object.");

  static Parameter<SimpleKTCut,Energy> interfaceMinKT
    ("MinKT",
     "The minimum allowed value of the transverse momentum of an outgoing "
     "parton.",
     &SimpleKTCut::theMinKT, GeV, 10.0*GeV, ZERO, Constants::MaxEnergy,
     true, false, Interface::limited,
     (ISFNK)0, (IGFNK)0, (IGFNK)0, &SimpleKTCut::maxKTMin, (IGFNK)0);


  static Parameter<SimpleKTCut,Energy> interfaceMaxKT
    ("MaxKT",
     "The maximum allowed value of the transverse momentum of an outgoing "
     "parton. Note that this cut does not increase the efficiency of the phase "
     "space generation, but is only applied as a post-cut.",
     &SimpleKTCut::theMaxKT, GeV, Constants::MaxEnergy, ZERO, ZERO,
     true, false, Interface::lowerlim,
     (ISFNK)0, (IGFNK)0,  &SimpleKTCut::minKTMax, (IGFNK)0, (IGFNK)0);


  static Parameter<SimpleKTCut,double> interfaceMinEta
    ("MinEta",
     "The minimum allowed pseudo-rapidity of an outgoing parton. "
     "The pseudo-rapidity is measured in the lab system.",
     &SimpleKTCut::theMinEta,
     -Constants::MaxRapidity, 0, Constants::MaxRapidity,
     true, false, Interface::upperlim,
     (ISFN)0, (IGFN)0, (IGFN)0, &SimpleKTCut::maxEtaMin, (IGFN)0);

  static Parameter<SimpleKTCut,double> interfaceMaxEta
    ("MaxEta",
     "The maximum allowed pseudo-rapidity of an outgoing parton. "
     "The pseudo-rapidity is measured in the lab system.",
     &SimpleKTCut::theMaxEta,
     Constants::MaxRapidity, -Constants::MaxRapidity, 0,
     true, false, Interface::lowerlim,
     (ISFN)0, (IGFN)0,  &SimpleKTCut::minEtaMax, (IGFN)0, (IGFN)0);

  static Reference<SimpleKTCut,MatcherBase> interfaceMatcher
    ("Matcher",
     "If non-null only particles matching this object will be affected "
     "by the cut.",
     &SimpleKTCut::theMatcher, true, false, true, true, false);

  interfaceMinKT.rank(10);
  interfaceMaxKT.rank(6);
  interfaceMinEta.rank(9);
  interfaceMaxEta.rank(8);
  interfaceMatcher.rank(7);
  interfaceMinKT.setHasDefault(false);
  interfaceMaxKT.setHasDefault(false);
  interfaceMinEta.setHasDefault(false);
  interfaceMaxEta.setHasDefault(false);

}

