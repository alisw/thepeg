// -*- C++ -*-
//
// MadGraphOneCut.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MadGraphOneCut class.
//

#include "MadGraphOneCut.h"
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

IBPtr MadGraphOneCut::clone() const {
  return new_ptr(*this);
}

IBPtr MadGraphOneCut::fullclone() const {
  return new_ptr(*this);
}

Energy MadGraphOneCut::minKT(tcPDPtr p) const {
  if ( cutType != PT || !checkType(p) ) return ZERO;
  return theCut*GeV;
}

double MadGraphOneCut::minEta(tcPDPtr p) const {
  if ( cutType != ETA || !checkType(p) ) return -Constants::MaxRapidity;
  return -theCut;
}

double MadGraphOneCut::maxEta(tcPDPtr p) const {
  if ( cutType != ETA || !checkType(p) ) return Constants::MaxRapidity;
  return theCut;
}

Energy MadGraphOneCut::minMaxKT(tcPDPtr p) const {
  if ( cutType != XPT || !checkType(p) ) return ZERO;
  return theCut*GeV;  
}

bool MadGraphOneCut::passCuts(tcCutsPtr parent,
			      tcPDPtr ptype, LorentzMomentum p) const {
  if ( !checkType(ptype) ) return true;
  if ( cutType == PT ) return p.perp() > theCut*GeV;
  if ( cutType == ETA ) {
    double y = p.rapidity() + parent->Y() + parent->currentYHat();
    return abs(p.mt()*sinh(y)) < p.perp()*sinh(theCut);
  }
  return true;
}

bool MadGraphOneCut::checkType(tcPDPtr p) const {
  switch ( abs(p->id()) ) {
  case ParticleID::d:
  case ParticleID::u:
  case ParticleID::s:
  case ParticleID::c:
  case ParticleID::g:
    return particleType == JET;
  case ParticleID::b:
    return particleType == JET || particleType == BOT;
  case ParticleID::gamma:
    return particleType == PHO;
  case ParticleID::eminus:
  case ParticleID::nu_e:
  case ParticleID::muminus:
  case ParticleID::nu_mu:
  case ParticleID::tauminus:
  case ParticleID::nu_tau:
    return particleType == LEP;
  default:
    return false;
  }
}

void MadGraphOneCut::persistentOutput(PersistentOStream & os) const {
  os << oenum(cutType) << oenum(particleType) << theCut;
}

void MadGraphOneCut::persistentInput(PersistentIStream & is, int) {
  is >> ienum(cutType) >> ienum(particleType) >> theCut;
}

ClassDescription<MadGraphOneCut> MadGraphOneCut::initMadGraphOneCut;
// Definition of the static class description member.

void MadGraphOneCut::Init() {

  static ClassDocumentation<MadGraphOneCut> documentation
    ("Objects of the MadGraphOneCut class can be created automatically by "
     "the MadGraphReader class when scanning event files for information "
     "about cuts. It is also possible to create objects by hand and use "
     "it as any other OneCutBase object.");

  static Switch<MadGraphOneCut,CutType> interfaceCutType
    ("CutType",
     "The type of cut this object will do.",
     &MadGraphOneCut::cutType, PT, true, false);
  static SwitchOption interfaceCutTypePT
    (interfaceCutType,
     "MinPT",
     "The minimum transverse momentum of a particle.",
     PT);
  static SwitchOption interfaceCutTypeMaxEta
    (interfaceCutType,
     "MaxEta",
     "The maximum (absolute value of) pseudo-rapidity of a particle.",
     ETA);
  static SwitchOption interfaceCutTypeMinMaxPT
    (interfaceCutType,
     "MinMaxPT",
     "The minimum transverse momentum of the particle with largest "
     "transverse momentum.",
     XPT);

  static Switch<MadGraphOneCut,PType> interfaceParticleType
    ("ParticleType",
     "The types of particles this cut is applied to.",
     &MadGraphOneCut::particleType, JET, true, false);
  static SwitchOption interfaceParticleTypeJets
    (interfaceParticleType,
     "Jets",
     "The cut applies only to coloured particles (jets).",
     JET);
  static SwitchOption interfaceParticleTypeLeptons
    (interfaceParticleType,
     "Leptons",
     "The cut applies only to leptons.",
     LEP);
  static SwitchOption interfaceParticleTypePhotons
    (interfaceParticleType,
     "Photons",
     "The cut applies only to photons.",
     PHO);
  static SwitchOption interfaceParticleTypeBottom
    (interfaceParticleType,
     "Bottom",
     "The cut applies only to bottom quarks.",
     BOT);

  static Parameter<MadGraphOneCut,double> interfaceCut
    ("Cut",
     "The value of the cut to be applied (in units of GeV in case of a "
     "transverse momentum).",
     &MadGraphOneCut::theCut, 0.0, 0.0, 0,
     true, false, Interface::lowerlim);

  interfaceCut.rank(10);
  interfaceCutType.rank(9);
  interfaceParticleType.rank(8);
  interfaceCut.setHasDefault(false);
  interfaceCutType.setHasDefault(false);
  interfaceParticleType.setHasDefault(false);

}

