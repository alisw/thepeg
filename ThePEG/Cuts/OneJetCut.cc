// -*- C++ -*-
//
// OneJetCut.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
// Copyright (C) 2009-2012 Simon Platzer
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneJetCut class.
//

#include "OneJetCut.h"
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

OneJetCut::OneJetCut()
  : ptMin(20.*GeV), yMin(-5), yMax(5) {}

OneJetCut::~OneJetCut() {}

void OneJetCut::describe() const {
  CurrentGenerator::log()
    << fullName() << " requesting one jet with:\n"
    << "pt(jet)/GeV > " << ptMin/GeV << " y(jet) > "
    << yMin << " y(jet) < " << yMax << "\n";
}

IBPtr OneJetCut::clone() const {
  return new_ptr(*this);
}

IBPtr OneJetCut::fullclone() const {
  return new_ptr(*this);
}      

void OneJetCut::persistentOutput(PersistentOStream & os) const {
  os << unresolvedMatcher << ounit(ptMin,GeV)
     << yMin << yMax;
}

void OneJetCut::persistentInput(PersistentIStream & is, int) {
  is >> unresolvedMatcher >> iunit(ptMin,GeV)
     >> yMin >> yMax;
}

bool OneJetCut::passCuts(tcCutsPtr, const tcPDVector & ptype,
			 const vector<LorentzMomentum> & p) const {
  tcPDVector::const_iterator ptit = ptype.begin();
  vector<LorentzMomentum>::const_iterator pit = p.begin();
  for ( ; ptit != ptype.end(); ++ptit, ++pit )
    if ( unresolvedMatcher->check(**ptit) ) {
      double y = pit->rapidity();
      if ( pit->perp() > ptMin &&
	   y > yMin && y < yMax ) {
	return true;
      }
    }
  return false;
}

ClassDescription<OneJetCut> OneJetCut::initOneJetCut;
// Definition of the static class description member.

void OneJetCut::Init() {

  static ClassDocumentation<OneJetCut> documentation
    ("OneJetsCut is a simple one-jet inclusive cut, requiring at least "
     "one jet above a certain pt in a given rapidity interval.");

  static Reference<OneJetCut,MatcherBase> interfaceUnresolvedMatcher
    ("UnresolvedMatcher",
     "A matcher identifying unresolved partons",
     &OneJetCut::unresolvedMatcher, false, false, true, false, false);


  static Parameter<OneJetCut,Energy> interfacePTMin
    ("PTMin",
     "Set the minimum pt required for the jet.",
     &OneJetCut::ptMin, GeV, 20.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);


  static Parameter<OneJetCut,double> interfaceYMin
    ("YMin",
     "Set the minimum rapidity required for the jet.",
     &OneJetCut::yMin, -5, 0, 0,
     false, false, Interface::nolimits);

  static Parameter<OneJetCut,double> interfaceYMax
    ("YMax",
     "Set the maximum rapidity required for the jet.",
     &OneJetCut::yMax, 5, 0, 0,
     false, false, Interface::nolimits);

}

