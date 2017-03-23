// -*- C++ -*-
//
// JetCuts.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
// Copyright (C) 2009-2012 Simon Platzer
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the JetCuts class.
//

#include "JetCuts.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace ThePEG;

JetCuts::JetCuts() 
  : theOrdering(orderPt) {}

JetCuts::~JetCuts() {}

void JetCuts::describe() const {
  CurrentGenerator::log() << name() << " cutting on jets ordered in "
			  << (ordering() == orderPt ? "pt " : "y ")
			  << " with regions:\n";
  for ( vector<Ptr<JetRegion>::ptr>::const_iterator r = jetRegions().begin();
	r != jetRegions().end(); ++r )
    (**r).describe();
  if ( !jetPairRegions().empty() ) {
    CurrentGenerator::log() << "and jet pair cuts:\n";
    for ( vector<Ptr<JetPairRegion>::ptr>::const_iterator r = jetPairRegions().begin();
	  r != jetPairRegions().end(); ++r )
      (**r).describe();
    for ( vector<Ptr<MultiJetRegion>::ptr>::const_iterator r = multiJetRegions().begin();
	  r != multiJetRegions().end(); ++r )
      (**r).describe();
  }
  if ( !jetVetoRegions().empty() ) {
    CurrentGenerator::log() << "vetoing jets inside:\n";
    for ( vector<Ptr<JetRegion>::ptr>::const_iterator r = jetVetoRegions().begin();
	  r != jetVetoRegions().end(); ++r )
      (**r).describe();
  }
  CurrentGenerator::log() << "\n";
}

IBPtr JetCuts::clone() const {
  return new_ptr(*this);
}

IBPtr JetCuts::fullclone() const {
  return new_ptr(*this);
}      

void JetCuts::persistentOutput(PersistentOStream & os) const {
  os << theUnresolvedMatcher 
     << theJetRegions << theJetVetoRegions 
     << theJetPairRegions << theMultiJetRegions 
     << theOrdering;
}

void JetCuts::persistentInput(PersistentIStream & is, int) {
  is >> theUnresolvedMatcher 
     >> theJetRegions >> theJetVetoRegions 
     >> theJetPairRegions >> theMultiJetRegions
     >> theOrdering;
}

struct PtLarger {
  inline bool operator()(const LorentzMomentum& a,
			 const LorentzMomentum& b) const {
    return a.perp() > b.perp();
  }
};

struct YLess {
  inline bool operator()(const LorentzMomentum& a,
			 const LorentzMomentum& b) const {
    return a.rapidity() < b.rapidity();
  }
};

bool JetCuts::passCuts(tcCutsPtr parent, const tcPDVector & ptype,
		       const vector<LorentzMomentum> & p) const {

  vector<LorentzMomentum> jets;
  tcPDVector::const_iterator ptypeit = ptype.begin();
  vector<LorentzMomentum>::const_iterator pit = p.begin();
  for ( ; ptypeit != ptype.end(); ++ptypeit, ++pit )
    if ( unresolvedMatcher()->check(**ptypeit) )
      jets.push_back(*pit);
  if ( ordering() == orderPt ) {
    sort(jets.begin(),jets.end(),PtLarger());
  } else if ( ordering() == orderY ) {
    sort(jets.begin(),jets.end(),YLess());
  } else assert(false);

  for ( vector<Ptr<JetRegion>::ptr>::const_iterator r = jetRegions().begin();
	r != jetRegions().end(); ++r )
    (**r).reset();

  set<int> matchedJets;

  for ( size_t k = 0; k < jets.size(); ++k ) {
    for ( vector<Ptr<JetRegion>::ptr>::const_iterator r = jetRegions().begin();
	  r != jetRegions().end(); ++r ) {
      if ( (**r).matches(parent,k+1,jets[k]) ) {
	matchedJets.insert(k+1);
	break;
      }
    }
  }

  double cutWeight = 1.0;
  bool pass = true;

  for ( vector<Ptr<JetRegion>::ptr>::const_iterator r = jetRegions().begin();
	r != jetRegions().end(); ++r ) {
    pass &= (**r).didMatch();
    cutWeight *= (**r).cutWeight();
    if ( !pass ) {
      parent->lastCutWeight(0.0);
      return false;
    }
  }

  for ( vector<Ptr<JetPairRegion>::ptr>::const_iterator r = jetPairRegions().begin();
	r != jetPairRegions().end(); ++r ) {
    pass &= (**r).matches(parent);
    cutWeight *= (**r).cutWeight();
    if ( !pass ) {
      parent->lastCutWeight(0.0);
      return false;
    }
  }

  for ( vector<Ptr<MultiJetRegion>::ptr>::const_iterator r = multiJetRegions().begin();
	r != multiJetRegions().end(); ++r ) {
    pass &= (**r).matches();
    cutWeight *= (**r).cutWeight();
    if ( !pass ) {
      parent->lastCutWeight(0.0);
      return false;
    }
  }

  if ( !jetVetoRegions().empty() ) {
    for ( size_t k = 0; k < jets.size(); ++k ) {
      if ( matchedJets.find(k+1) != matchedJets.end() )
	continue;
      for ( vector<Ptr<JetRegion>::ptr>::const_iterator r = jetVetoRegions().begin();
	    r != jetVetoRegions().end(); ++r ) {
	pass &= !(**r).matches(parent,k+1,jets[k]);
	cutWeight *= 1. - (**r).cutWeight();
	if ( !pass ) {
	  parent->lastCutWeight(0.0);
	  return false;
	}
      }
    }
  }

  parent->lastCutWeight(cutWeight);

  return pass;

}

ClassDescription<JetCuts> JetCuts::initJetCuts;
// Definition of the static class description member.

void JetCuts::Init() {

  static ClassDocumentation<JetCuts> documentation
    ("JetCuts combines various JetRegion and JetPairRegion objects into a "
     "cut object.");

  static Reference<JetCuts,MatcherBase> interfaceUnresolvedMatcher
    ("UnresolvedMatcher",
     "A matcher identifying unresolved partons",
     &JetCuts::theUnresolvedMatcher, false, false, true, false, false);

  static RefVector<JetCuts,JetRegion> interfaceJetRegions
    ("JetRegions",
     "The jet regions to be used.",
     &JetCuts::theJetRegions, -1, false, false, true, false, false);

  static RefVector<JetCuts,JetRegion> interfaceJetVetoRegions
    ("JetVetoRegions",
     "The jet veto regions to be used.",
     &JetCuts::theJetVetoRegions, -1, false, false, true, false, false);

  static RefVector<JetCuts,JetPairRegion> interfaceJetPairRegions
    ("JetPairRegions",
     "The jet pair regions to be used.",
     &JetCuts::theJetPairRegions, -1, false, false, true, false, false);

  static RefVector<JetCuts,MultiJetRegion> interfaceMultiJetRegions
    ("MultiJetRegions",
     "The multi jet regions to be used.",
     &JetCuts::theMultiJetRegions, -1, false, false, true, false, false);


  static Switch<JetCuts,int> interfaceOrdering
    ("Ordering",
     "The ordering to apply on jets.",
     &JetCuts::theOrdering, orderPt, false, false);
  static SwitchOption interfaceOrderingOrderPt
    (interfaceOrdering,
     "OrderPt",
     "Order in decreasing transverse momentum.",
     orderPt);
  static SwitchOption interfaceOrderingOrderY
    (interfaceOrdering,
     "OrderY",
     "Order in rapidity.",
     orderY);

}

