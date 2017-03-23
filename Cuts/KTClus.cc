// -*- C++ -*-
//
// KTClus.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KTClus class.
//

#include "KTClus.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace ThePEG;

void KTClus::describe() const {
  CurrentGenerator::log() 
    << fullName() << ":\n"
    << "Cut = " << theCut/GeV << " GeV\n\n";
}

IBPtr KTClus::clone() const {
  return new_ptr(*this);
}

IBPtr KTClus::fullclone() const {
  return new_ptr(*this);
}

Energy KTClus::minKTClus(tcPDPtr pi, tcPDPtr pj) const {
  if ( onlyJets &&
       ( ( pi && !pi->coloured() ) ||
	 ( pj && !pj->coloured() ) ) ) return ZERO;
  if ( !pi && !pj ) return ZERO;
  return theCut;
}

Energy2 KTClus::minSij(tcPDPtr, tcPDPtr) const {
  return ZERO;
}

Energy2 KTClus::minTij(tcPDPtr, tcPDPtr) const {
  return ZERO;
}

double KTClus::minDeltaR(tcPDPtr, tcPDPtr) const {
  return 0.0;
}

double KTClus::minDurham(tcPDPtr, tcPDPtr) const {
  return 0.0;
}

bool KTClus::passCuts(tcCutsPtr, tcPDPtr pitype, tcPDPtr pjtype,
		      LorentzMomentum pi, LorentzMomentum pj,
		      bool inci, bool incj) const {
  if ( onlyJets &&
       ( ( !inci && pitype && !pitype->coloured() ) ||
	 ( !incj && pjtype && !pjtype->coloured() ) ) ) return true;
  if ( inci && incj ) return true;
  else if ( inci ) {
    if ( pj.perp() <= theCut ) return false;
  }
  else if ( incj ) {
    if ( pi.perp() <= theCut ) return false;
  }
  else {
    double deta2 = sqr(pi.eta() - pj.eta());
    double dphi = abs(pi.phi() - pj.phi());
    if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
    double dr = sqrt(deta2 + sqr(dphi));
    if ( min(pi.perp(), pj.perp())*dr <= theCut ) return false;
  }
  return true;
}

void KTClus::persistentOutput(PersistentOStream & os) const {
  os << ounit(theCut, GeV) << onlyJets;
}

void KTClus::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theCut, GeV) >> onlyJets;
}

ClassDescription<KTClus> KTClus::initKTClus;
// Definition of the static class description member.

void KTClus::Init() {

  static ClassDocumentation<KTClus> documentation
    ("This clas implements the cuts relevant for the "
     "\\f$\\Delta R\\f$-measure in the longitudinally invariant "
     "kt-algorithm. By default the cut is only applied to coloured "
     "particles, but optionally it may be applied to all particle types.");

  static Parameter<KTClus,Energy> interfaceCut
    ("Cut",
     "The minimum allowed value of the longitudinally invariant "
     "\\f$k_\\perp\\f$-algorithms distance measure. This is defined as "
     "\\f$\\min(p_{\\perp i}, p_{\\perp j})"
     "\\sqrt{\\Delta\\eta_{ij}^2+\\Delta\\phi_{ij}^2}\\f$ for two outgoing "
     "partons, or simply \\f$p_{\\perp i}\\f$ or \\f$p_{\\perp j}\\f$ for a "
     "single outgoing parton.",
     &KTClus::theCut, GeV, 10.0*GeV, ZERO, ZERO,
     true, false, Interface::lowerlim);
  interfaceCut.setHasDefault(false);

  static Switch<KTClus,bool> interfaceOnlyJets
    ("OnlyJets",
     "If true, only apply the cut to pairs of coloured particles.",
     &KTClus::onlyJets, true, true, false);
  static SwitchOption interfaceOnlyJetsOnlyJets
    (interfaceOnlyJets,
     "OnlyJets",
     "The cut is only applied to pairs of coloured particles.",
     true);
  static SwitchOption interfaceOnlyJetsAll
    (interfaceOnlyJets,
     "All",
     "The cut is applied to all particle types.",
     false);
  interfaceOnlyJets.setHasDefault(false);

  interfaceCut.rank(10);
  interfaceOnlyJets.rank(9);

}

