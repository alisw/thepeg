// -*- C++ -*-
//
// NJetsCut.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
// Copyright (C) 2009-2012 Simon Platzer
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NJetsCut class.
//

#include "NJetsCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace ThePEG;

NJetsCut::NJetsCut()
  : nJetsMin(0), nJetsMax(-1) {}

NJetsCut::~NJetsCut() {}

void NJetsCut::describe() const {
  CurrentGenerator::log()
    << fullName() << ": requires ";
  if ( nJetsMin > 0 )
    CurrentGenerator::log() << "at least "
			    << nJetsMin;
  if ( nJetsMax > 0 )
    CurrentGenerator::log() << " and at most "
			    << nJetsMax;
  CurrentGenerator::log() << " jets.\n";
}

IBPtr NJetsCut::clone() const {
  return new_ptr(*this);
}

IBPtr NJetsCut::fullclone() const {
  return new_ptr(*this);
}      

void NJetsCut::persistentOutput(PersistentOStream & os) const {
  os << unresolvedMatcher
     << nJetsMin << nJetsMax;
}

void NJetsCut::persistentInput(PersistentIStream & is, int) {
  is >> unresolvedMatcher
     >> nJetsMin >> nJetsMax;
}

bool NJetsCut::passCuts(tcCutsPtr, const tcPDVector & ptype,
			 const vector<LorentzMomentum> & p) const {
  tcPDVector::const_iterator ptit = ptype.begin();
  vector<LorentzMomentum>::const_iterator pit = p.begin();
  int njets = 0;
  for ( ; ptit != ptype.end(); ++ptit, ++pit )
    if ( unresolvedMatcher->check(**ptit) ) {
      ++njets;
    }
  if ( nJetsMax > 0 )
    return njets >= nJetsMin && njets <= nJetsMax;
  return njets >= nJetsMin;
}

ClassDescription<NJetsCut> NJetsCut::initNJetsCut;
// Definition of the static class description member.

void NJetsCut::Init() {

  static ClassDocumentation<NJetsCut> documentation
    ("NJetsCut is a simple cut on jet multiplicity.");

  static Reference<NJetsCut,MatcherBase> interfaceUnresolvedMatcher
    ("UnresolvedMatcher",
     "A matcher identifying unresolved partons",
     &NJetsCut::unresolvedMatcher, false, false, true, false, false);


  static Parameter<NJetsCut,int> interfaceNJetsMin
    ("NJetsMin",
     "The minimum number of jets required.",
     &NJetsCut::nJetsMin, 0, 0, 0,
     false, false, Interface::lowerlim);


  static Parameter<NJetsCut,int> interfaceNJetsMax
    ("NJetsMax",
     "The maximum number of jets allowed. If -1 no limit is imposed.",
     &NJetsCut::nJetsMax, -1, -1, 0,
     false, false, Interface::lowerlim);

}

