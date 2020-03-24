// -*- C++ -*-
//
// MadGraphTwoCut.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
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
  if ( !checkType(pi, pj) || cutType != Cut::INVMASS ) return ZERO;
  return sqr(theCut*GeV);
}

Energy2 MadGraphTwoCut::minTij(tcPDPtr, tcPDPtr) const {
  return ZERO;
}

double MadGraphTwoCut::minDeltaR(tcPDPtr pi, tcPDPtr pj) const {
  if ( !checkType(pi, pj) || cutType != Cut::DELTAR ) return 0.0;
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
  if ( cutType == Cut::INVMASS ) return (pi + pj).m2() > sqr(theCut*GeV);
  if ( cutType == Cut::DELTAR ) {
    double deta2 = sqr(pi.eta() - pj.eta());
    double dphi = abs(pi.phi() - pj.phi());
    if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
    return sqrt(deta2 + sqr(dphi)) > theCut;
  }
  return true;
}

MadGraphTwoCut::P MadGraphTwoCut::getType(tcPDPtr p) const {
  switch ( abs(p->id()) ) {
  case ParticleID::d:
  case ParticleID::u:
  case ParticleID::s:
  case ParticleID::c:
  case ParticleID::g:
    return P::JET;
  case ParticleID::b:
    return P::BOT;
  case ParticleID::gamma:
    return P::PHO;
  case ParticleID::eminus:
  case ParticleID::nu_e:
  case ParticleID::muminus:
  case ParticleID::nu_mu:
  case ParticleID::tauminus:
  case ParticleID::nu_tau:
    return P::LEP;
  default:
    return P::NOT;
  }
}

bool MadGraphTwoCut::checkType(tcPDPtr pi, tcPDPtr pj) const {
  switch ( pairType ) {
  case PP::JETJET:
    return getType(pi) == P::JET && getType(pj) == P::JET;
  case PP::LEPLEP:
    if ( getType(pi) != P::LEP || getType(pj) != P::LEP ) return false;
    if ( cutType == Cut::DELTAR ) return true;
    // Special treatment for Cut::INVMASS.
    if ( pi->id()*pj->id() >= 0 ) return false;
    // OK we have a lepton-anti-lepton pair. I it the same lepton
    if ( pi->id() == -pj->id() ) return true;
    // NO, well is it the same family?
    if ( max(abs(pi->id()), abs(pj->id()))%2 ) return false;
    return abs(pi->id() + pj->id()) == 1 ;
  case PP::PHOPHO:
    return getType(pi) == P::PHO && getType(pj) == P::PHO;
  case PP::BOTBOT:
    return getType(pi) == P::BOT && getType(pj) == P::BOT;
  case PP::BOTJET:
    return ( getType(pi) == P::BOT && getType(pj) == P::JET ) ||
           ( getType(pi) == P::JET && getType(pj) == P::BOT );
  case PP::PHOJET:
    return ( getType(pi) == P::PHO && getType(pj) == P::JET ) ||
           ( getType(pi) == P::JET && getType(pj) == P::PHO );
  case PP::JETLEP:
    return ( getType(pi) == P::LEP && getType(pj) == P::JET ) ||
           ( getType(pi) == P::JET && getType(pj) == P::LEP );
  case PP::PHOBOT:
    return ( getType(pi) == P::PHO && getType(pj) == P::BOT ) ||
           ( getType(pi) == P::BOT && getType(pj) == P::PHO );
  case PP::BOTLEP:
    return ( getType(pi) == P::BOT && getType(pj) == P::LEP ) ||
           ( getType(pi) == P::LEP && getType(pj) == P::BOT );
  case PP::PHOLEP:
    return ( getType(pi) == P::PHO && getType(pj) == P::LEP ) ||
           ( getType(pi) == P::LEP && getType(pj) == P::PHO );
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

  static Switch<MadGraphTwoCut,Cut> interfaceCutType
    ("CutType",
     "The kind of cut this object will do.",
     &MadGraphTwoCut::cutType, Cut::DELTAR, true, false);
  static SwitchOption interfaceCutTypeInvariantMass
    (interfaceCutType,
     "InvariantMass",
     "The minimum invariant mass of two particles.",
     Cut::INVMASS);
  static SwitchOption interfaceCutTypeDeltaR
    (interfaceCutType,
     "DeltaR",
     "The minimum pseudo-rapidity--azimuth-angle distance between two "
     "particles.",
     Cut::DELTAR);

  static Switch<MadGraphTwoCut,PP> interfacePairType
    ("PairType",
     "The type of particle pairs this cut is applied to.",
     &MadGraphTwoCut::pairType, PP::JETJET, true, false);
  static SwitchOption interfacePairTypeJetJet
    (interfacePairType,
     "JetJet",
     "The cut applies only to pairs of coloured particles (jets).",
     PP::JETJET);
  static SwitchOption interfacePairTypeLeptonLepton
    (interfacePairType,
     "LeptonLepton",
     "The cut applies only to lepton pairs (in case of invariant mass, "
     "lepton--anti-lepton pairs of same flavour).",
     PP::LEPLEP);
  static SwitchOption interfacePairTypePhotonPhoton
    (interfacePairType,
     "PhotonPhoton",
     "The cut applies only to pairs photons.",
     PP::PHOPHO);
  static SwitchOption interfacePairTypeBottomPairs
    (interfacePairType,
     "BottomPairs",
     "The cut applies only to pairs of bottom quarks.",
     PP::BOTBOT);
  static SwitchOption interfacePairTypeJetBottom
    (interfacePairType,
     "JetBottom",
     "The cut applies only to bottom quarks paired with another coloured "
     "particle (jet).",
     PP::BOTJET);
  static SwitchOption interfacePairTypePhotonJet
    (interfacePairType,
     "PhotonJet",
     "The cut applies only to a photon paired with a coloured particle (jet).",
     PP::PHOJET);
  static SwitchOption interfacePairTypeJetLepton
    (interfacePairType,
     "JetLepton",
     "The cut applies only to a coloured particle (jet) paired with a lepton.",
     PP::JETLEP);
  static SwitchOption interfacePairTypePhotonBottom
    (interfacePairType,
     "PhotonBottom",
     "The cut applies only to a photon paired with a bottom quark.",
     PP::PHOBOT);
  static SwitchOption interfacePairTypeBottomLepton
    (interfacePairType,
     "BottomLepton",
     "The cut applies only to bottom quarks paired with a lepton.",
     PP::BOTLEP);
  static SwitchOption interfacePairTypePhotonLepton
    (interfacePairType,
     "PhotonLepton",
     "The cut applies only to a photon paired with a lepton.",
     PP::PHOLEP);

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

