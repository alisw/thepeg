// -*- C++ -*-
//
// MadGraphTwoCut.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MadGraphTwoCut class.
//

#include "MadGraphTwoCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

IBPtr MadGraphTwoCut::clone() const {
  return new_ptr(*this);
}

IBPtr MadGraphTwoCut::fullclone() const {
  return new_ptr(*this);
}

Energy2 MadGraphTwoCut::minSij(tcPDPtr pi, tcPDPtr pj) const {
  if ( !checkType(pi, pj) || cutType != INVMASS ) return ZERO;
  return sqr(theCut*GeV);
}

Energy2 MadGraphTwoCut::minTij(tcPDPtr, tcPDPtr) const {
  return ZERO;
}

double MadGraphTwoCut::minDeltaR(tcPDPtr pi, tcPDPtr pj) const {
  if ( !checkType(pi, pj) || cutType != DELTAR ) return 0.0;
  return theCut;
}

Energy MadGraphTwoCut::minKTClus(tcPDPtr, tcPDPtr) const {
  return ZERO;
}

double MadGraphTwoCut::minDurham(tcPDPtr, tcPDPtr) const {
  return 0.0;
}

bool MadGraphTwoCut::passCuts(tcCutsPtr, tcPDPtr pitype, tcPDPtr pjtype,
			      LorentzMomentum pi, LorentzMomentum pj,
			      bool inci, bool incj) const {
  if ( inci || incj || !checkType(pitype, pjtype) ) return true;
  if ( cutType == INVMASS ) return (pi + pj).m2() > sqr(theCut*GeV);
  if ( cutType == DELTAR ) {
    double deta2 = sqr(pi.eta() - pj.eta());
    double dphi = abs(pi.phi() - pj.phi());
    if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
    return sqrt(deta2 + sqr(dphi)) > theCut;
  }
  return true;
}

MadGraphTwoCut::PType MadGraphTwoCut::getType(tcPDPtr p) const {
  switch ( abs(p->id()) ) {
  case ParticleID::d:
  case ParticleID::u:
  case ParticleID::s:
  case ParticleID::c:
  case ParticleID::g:
    return JET;
  case ParticleID::b:
    return BOT;
  case ParticleID::gamma:
    return PHO;
  case ParticleID::eminus:
  case ParticleID::nu_e:
  case ParticleID::muminus:
  case ParticleID::nu_mu:
  case ParticleID::tauminus:
  case ParticleID::nu_tau:
    return LEP;
  default:
    return NOT;
  }
}

bool MadGraphTwoCut::checkType(tcPDPtr pi, tcPDPtr pj) const {
  switch ( pairType ) {
  case JETJET:
    return getType(pi) == JET && getType(pj) == JET;
  case LEPLEP:
    if ( getType(pi) != LEP || getType(pj) != LEP ) return false;
    if ( cutType == DELTAR ) return true;
    // Special treatment for INVMASS.
    if ( pi->id()*pj->id() >= 0 ) return false;
    // OK we have a lepton-anti-lepton pair. I it the same lepton
    if ( pi->id() == -pj->id() ) return true;
    // NO, well is it the same family?
    if ( max(abs(pi->id()), abs(pj->id()))%2 ) return false;
    return abs(pi->id() + pj->id()) == 1 ;
  case PHOPHO:
    return getType(pi) == PHO && getType(pj) == PHO;
  case BOTBOT:
    return getType(pi) == BOT && getType(pj) == BOT;
  case BOTJET:
    return ( getType(pi) == BOT && getType(pj) == JET ) ||
           ( getType(pi) == JET && getType(pj) == BOT );
  case PHOJET:
    return ( getType(pi) == PHO && getType(pj) == JET ) ||
           ( getType(pi) == JET && getType(pj) == PHO );
  case JETLEP:
    return ( getType(pi) == LEP && getType(pj) == JET ) ||
           ( getType(pi) == JET && getType(pj) == LEP );
  case PHOBOT:
    return ( getType(pi) == PHO && getType(pj) == BOT ) ||
           ( getType(pi) == BOT && getType(pj) == PHO );
  case BOTLEP:
    return ( getType(pi) == BOT && getType(pj) == LEP ) ||
           ( getType(pi) == LEP && getType(pj) == BOT );
  case PHOLEP:
    return ( getType(pi) == PHO && getType(pj) == LEP ) ||
           ( getType(pi) == LEP && getType(pj) == PHO );
  }
  return false;
}

void MadGraphTwoCut::persistentOutput(PersistentOStream & os) const {
  os << oenum(cutType) << oenum(pairType) << theCut;
}

void MadGraphTwoCut::persistentInput(PersistentIStream & is, int) {
  is >> ienum(cutType) >> ienum(pairType) >> theCut;
}

ClassDescription<MadGraphTwoCut> MadGraphTwoCut::initMadGraphTwoCut;
// Definition of the static class description member.

void MadGraphTwoCut::Init() {

  static ClassDocumentation<MadGraphTwoCut> documentation
    ("Objects of the MadGraphTwoCut class can be created automatically by "
     "the MadGraphReader class when scanning event files for information "
     "about cuts. It is also possible to create objects by hand and use "
     "it as any other MadGraphTwoCut object.");

  static Switch<MadGraphTwoCut,CutType> interfaceCutType
    ("CutType",
     "The kind of cut this object will do.",
     &MadGraphTwoCut::cutType, DELTAR, true, false);
  static SwitchOption interfaceCutTypeInvariantMass
    (interfaceCutType,
     "InvariantMass",
     "The minimum invariant mass of two particles.",
     INVMASS);
  static SwitchOption interfaceCutTypeDeltaR
    (interfaceCutType,
     "DeltaR",
     "The minimum pseudo-rapidity--azimuth-angle distance between two "
     "particles.",
     DELTAR);

  static Switch<MadGraphTwoCut,PPType> interfacePairType
    ("PairType",
     "The type of particle pairs this cut is applied to.",
     &MadGraphTwoCut::pairType, JETJET, true, false);
  static SwitchOption interfacePairTypeJetJet
    (interfacePairType,
     "JetJet",
     "The cut applies only to pairs of coloured particles (jets).",
     JETJET);
  static SwitchOption interfacePairTypeLeptonLepton
    (interfacePairType,
     "LeptonLepton",
     "The cut applies only to lepton pairs (in case of invariant mass, "
     "lepton--anti-lepton pairs of same flavour).",
     LEPLEP);
  static SwitchOption interfacePairTypePhotonPhoton
    (interfacePairType,
     "PhotonPhoton",
     "The cut applies only to pairs photons.",
     PHOPHO);
  static SwitchOption interfacePairTypeBottomPairs
    (interfacePairType,
     "BottomPairs",
     "The cut applies only to pairs of bottom quarks.",
     BOTBOT);
  static SwitchOption interfacePairTypeJetBottom
    (interfacePairType,
     "JetBottom",
     "The cut applies only to bottom quarks paired with another coloured "
     "particle (jet).",
     BOTJET);
  static SwitchOption interfacePairTypePhotonJet
    (interfacePairType,
     "PhotonJet",
     "The cut applies only to a photon paired with a coloured particle (jet).",
     PHOJET);
  static SwitchOption interfacePairTypeJetLepton
    (interfacePairType,
     "JetLepton",
     "The cut applies only to a coloured particle (jet) paired with a lepton.",
     JETLEP);
  static SwitchOption interfacePairTypePhotonBottom
    (interfacePairType,
     "PhotonBottom",
     "The cut applies only to a photon paired with a bottom quark.",
     PHOBOT);
  static SwitchOption interfacePairTypeBottomLepton
    (interfacePairType,
     "BottomLepton",
     "The cut applies only to bottom quarks paired with a lepton.",
     BOTLEP);
  static SwitchOption interfacePairTypePhotonLepton
    (interfacePairType,
     "PhotonLepton",
     "The cut applies only to a photon paired with a lepton.",
     PHOLEP);

  static Parameter<MadGraphTwoCut,double> interfaceCut
    ("Cut",
     "The value of the cut to be applied (in units of GeV in case of "
     "minimum invariant mass).",
     &MadGraphTwoCut::theCut, 0.0, 0.0, 0,
     true, false, Interface::lowerlim);

  interfaceCut.rank(10);
  interfaceCutType.rank(9);
  interfacePairType.rank(8);
  interfaceCut.setHasDefault(false);
  interfaceCutType.setHasDefault(false);
  interfacePairType.setHasDefault(false);

}

