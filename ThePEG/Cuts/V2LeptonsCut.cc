// -*- C++ -*-
//
// V2LeptonsCut.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the V2LeptonsCut class.
//

#include "V2LeptonsCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace ThePEG;

V2LeptonsCut::~V2LeptonsCut() {}

void V2LeptonsCut::describe() const {
  CurrentGenerator::log() 
    << fullName() << ":\n"
    << "M = " << theMinM/GeV << " .. " << theMaxM/GeV << " GeV\n\n";
}

IBPtr V2LeptonsCut::clone() const {
  return new_ptr(*this);
}

IBPtr V2LeptonsCut::fullclone() const {
  return new_ptr(*this);
}

Energy2 V2LeptonsCut::minS(const tcPDVector & pv) const {
  if ( pv.size() != 2 ) return ZERO;
  if ( !checkTypes(pv[0]->id(), pv[1]->id()) ) return ZERO;
  return sqr(theMinM);
}

Energy2 V2LeptonsCut::maxS(const tcPDVector & pv) const {
  if ( pv.size() != 2 ) return Constants::MaxEnergy2;
  if ( !checkTypes(pv[0]->id(), pv[1]->id()) ) return Constants::MaxEnergy2;
  return sqr(theMaxM);
}

bool V2LeptonsCut::passCuts(tcCutsPtr, const tcPDVector & ptype,
	      const vector<LorentzMomentum> & p) const {
  for ( int i = 0, N = ptype.size() - 1; i < N; ++i )
    for ( int j = i + 1, M = ptype.size(); j < M; ++j ) {
      if ( !checkTypes(ptype[i]->id(), ptype[j]->id()) ) continue;
      Energy2 s = (p[i] + p[j]).m2();
      if ( s <= sqr(theMinM) || s >= sqr(theMaxM) ) return false;
    }
  return true;
}
      

int V2LeptonsCut::family(long id) const {
  switch ( abs(id) ) {
  case ParticleID::eminus:
  case ParticleID::nu_e:
    return electron;
  case ParticleID::muminus:
  case ParticleID::nu_mu:
    return muon;
  case ParticleID::tauminus:
  case ParticleID::nu_tau:
    return tau;
  }
  return 0;
}

bool V2LeptonsCut::checkTypes(long id1, long id2) const {
  // Must be particle anti-particle pair;
  if ( id1*id2 >= 0 ) return false;

  // Check that we have leptons, the families are the same and matches
  // the chosen ones.
  int fam1 = family(id1);
  if ( !fam1 ) return false;
  int fam2 = family(id2);
  if ( fam2 != fam1 || !(theFamilies&fam1) ) return false;

  // Check charge combination.
  int ccomb;
  if ( (id1%2) && (id2%2) ) ccomb = posneg;
  else if ( id1%2 ) {
    if ( id1 > 0 ) ccomb = negneu;
    else ccomb = posneu;
  }
  else if ( id2%2 ) {
    if ( id2 > 0 ) ccomb = negneu;
    else ccomb = posneu;
  }
  else
    ccomb = neuneu;

  return (theCComb&ccomb);
}

void V2LeptonsCut::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMinM, GeV) << ounit(theMaxM, GeV) << theFamilies << theCComb;
}

void V2LeptonsCut::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMinM, GeV) >> iunit(theMaxM, GeV) >> theFamilies >> theCComb;
}

Energy V2LeptonsCut::maxMinM() const {
  return theMaxM;
}

Energy V2LeptonsCut::minMaxM() const {
  return theMinM;
}

ClassDescription<V2LeptonsCut> V2LeptonsCut::initV2LeptonsCut;
// Definition of the static class description member.

void V2LeptonsCut::Init() {

  static ClassDocumentation<V2LeptonsCut> documentation
    ("This class inherits from MultiCutBase and describes cuts on the "
     "invariant mass of two final state leptons corresponding to the decay "
     "of a vector boson. It can be used when generating matrix elements to "
     "avoid the long tails of the resonance.");

  static Parameter<V2LeptonsCut,Energy> interfaceMinM
    ("MinM",
     "The minimum allowed invariant mass of the matched lepton pair.",
     &V2LeptonsCut::theMinM, GeV, 70.0*GeV, ZERO, Constants::MaxEnergy,
     true, false, Interface::limited,
     0, 0, 0, &V2LeptonsCut::maxMinM, 0);

  static Parameter<V2LeptonsCut,Energy> interfaceMaxM
    ("MaxM",
     "The maximum allowed invariant mass of the matched lepton pair.",
     &V2LeptonsCut::theMaxM, GeV, 90.0*GeV, ZERO, ZERO,
     true, false, Interface::lowerlim,
     0, 0, &V2LeptonsCut::minMaxM, 0, 0);

  static Switch<V2LeptonsCut,int> interfaceFamilies
    ("Families",
     "The different lepton families for which this cut should apply.",
     &V2LeptonsCut::theFamilies, electron|muon, true, false);
  static SwitchOption interfaceFamiliesElectron
    (interfaceFamilies,
     "Electron",
     "Only apply cut to electrons and electron neutrinos.",
     electron);
  static SwitchOption interfaceFamiliesMuon
    (interfaceFamilies,
     "Muon",
     "Only apply cut to muons and muon neutrinos.",
     muon);
  static SwitchOption interfaceFamiliesTau
    (interfaceFamilies,
     "Tau",
     "Only apply cut to taus and tau neutrinos.",
     tau);
  static SwitchOption interfaceFamiliesElectronMuon
    (interfaceFamilies,
     "ElectronMuon",
     "Only apply cut to electron and muon leptons.",
     electron|muon);
  static SwitchOption interfaceFamiliesAll
    (interfaceFamilies,
     "All",
     "Apply cut to all lepton families.",
     electron|muon|tau);

  static Switch<V2LeptonsCut,int> interfaceCComb
    ("CComb",
     "The charge combination of the lepton pair on which to cut.",
     &V2LeptonsCut::theCComb, posneu|negneu, true, false);
  static SwitchOption interfaceCCombAll
    (interfaceCComb,
     "All",
     "Cut on all relevant charge combinations.",
     posneg|negneu|posneu|neuneu);
  static SwitchOption interfaceCCombWplus
    (interfaceCComb,
     "Wplus",
     "Cut on positive lepton neutrin pairs.",
     posneu);
  static SwitchOption interfaceCCombWminus
    (interfaceCComb,
     "Wminus",
     "Cut on negative lepton anti-neutrin pairs.",
     negneu);
  static SwitchOption interfaceCCombW
    (interfaceCComb,
     "W",
     "Cut on charged lepton neutrino pairs.",
     posneu|negneu);
  static SwitchOption interfaceCCombGamma
    (interfaceCComb,
     "Gamma",
     "Cut on charged lepton anti-lepton pairs.",
     posneg);
  static SwitchOption interfaceCCombZ
    (interfaceCComb,
     "Z",
     "Cut on lepton anti-lepton pairs.",
     neuneu|posneg);
  static SwitchOption interfaceCCombZneutrinos
    (interfaceCComb,
     "Zneutrinos",
     "Cut on neutrino anti-neutrino pairs.",
     neuneu);

  interfaceMinM.rank(10);
  interfaceMaxM.rank(9);
  interfaceMinM.setHasDefault(false);
  interfaceMaxM.setHasDefault(false);
  interfaceCComb.setHasDefault(false);
  interfaceFamilies.setHasDefault(false);
}

