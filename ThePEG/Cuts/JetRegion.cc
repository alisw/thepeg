// -*- C++ -*-
//
// JetRegion.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2007 Leif Lonnblad
// Copyright (C) 2009-2012 Simon Platzer
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the JetRegion class.
//

#include "JetRegion.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Repository/CurrentGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

JetRegion::JetRegion()
  : thePtMin(0.*GeV), thePtMax(Constants::MaxEnergy),
    theDidMatch(false), theLastNumber(0), 
    theFuzzy(false), theCutWeight(1.0),
    theEnergyCutWidth(1.0*GeV), theRapidityCutWidth(0.1) {}

JetRegion::~JetRegion() {}

IBPtr JetRegion::clone() const {
  return new_ptr(*this);
}

IBPtr JetRegion::fullclone() const {
  return new_ptr(*this);
}

string JetRegion::doYRange(string in) {
  istringstream ins(StringUtils::stripws(in));
  theYRanges.clear();
  while ( !ins.eof() ) {
    double first, second;
    ins >> first;
    if ( ins.eof() )
      throw Exception() << "need an even number n of values to define n/2 rapidity intervals"
			<< Exception::abortnow;
    ins >> second;
    if ( first > second )
      swap(first,second);
    theYRanges.push_back(make_pair(first,second));
  }
  return "";
}

void JetRegion::describe() const {

  CurrentGenerator::log()
    << "JetRegion '" << name() << "' matching ";
  if ( accepts().empty() )
    CurrentGenerator::log() << "any jets ";
  else {
    CurrentGenerator::log() << "jets ";
    for ( vector<int>::const_iterator k = accepts().begin();
	  k != accepts().end(); ++k ) {
      CurrentGenerator::log() << "#" << *k;
      if ( k != --accepts().end() )
	CurrentGenerator::log() << ", ";
      else
	CurrentGenerator::log() << " ";
    }
  }
  CurrentGenerator::log() << " within:\n";

  CurrentGenerator::log() 
    << "pt  = " << ptMin()/GeV << " .. " << ptMax()/GeV << " GeV\n";

  for ( vector<pair<double,double> >::const_iterator r = yRanges().begin();
	r != yRanges().end(); ++r ) {
    CurrentGenerator::log() << "y   = " << r->first << " .. " << r->second << "\n";
  }

}

double step(double r) {
  if ( r < -0.5 )
    return 0.0;
  if ( r > 0.5 )
    return 1.0;
  return r+0.5;
}

bool JetRegion::lessThanEnergy(Energy a, Energy b, double& weight) const {
  if ( !fuzzy() ) {
    if ( a < b ) {
      weight = 1.0;
      return true;
    }
    weight = 0.0;
    return false;
  }
  double w = step((b-a)/theEnergyCutWidth);
  if ( w == 0.0 ) {
    weight = 0.0;
    return false;
  }
  weight *= w;
  return true;
}

bool JetRegion::lessThanRapidity(double a, double b, double& weight) const {
  if ( !fuzzy() ) {
    if ( a < b ) {
      weight = 1.0;
      return true;
    }
    weight = 0.0;
    return false;
  }
  double w = step((b-a)/theRapidityCutWidth);
  if ( w == 0.0 ) {
    weight = 0.0;
    return false;
  }
  weight *= w;
  return true;
}

bool JetRegion::matches(tcCutsPtr parent, int n, const LorentzMomentum& p,
			double yHat) {

  // one jet region can only contain one jet
  if ( didMatch() )
    return false;

  if ( !accepts().empty() && find(accepts().begin(),accepts().end(),n) == accepts().end() )
    return false;

  theCutWeight = 1.0;

  if ( !(lessThanEnergy(ptMin(),p.perp(),theCutWeight) &&
	 lessThanEnergy(p.perp(),ptMax(),theCutWeight)) ) {
    theCutWeight = 0.0;
    theDidMatch = false;
    return false;
  }

  double currentY = parent ? parent->currentYHat() : yHat;

  bool inRange = false || yRanges().empty();
  for ( vector<pair<double,double> >::const_iterator r = yRanges().begin();
	r != yRanges().end(); ++r ) {
    double rangeWeight = 1.0;
    if ( lessThanRapidity(r->first,p.rapidity() + currentY,rangeWeight) &&
	 lessThanRapidity(p.rapidity() + currentY,r->second,rangeWeight) ) {
      theCutWeight *= rangeWeight;
      inRange = true;
      break;
    }
  }
  if ( !inRange ) {
    theCutWeight = 0.0;
    theDidMatch = false;
    return false;
  }

  theDidMatch = true;
  theLastNumber = n;
  theLastMomentum = p;

  return true;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void JetRegion::persistentOutput(PersistentOStream & os) const {
  os << ounit(thePtMin,GeV) << ounit(thePtMax,GeV)
     << theYRanges << theAccepts << theFuzzy << theCutWeight
     << ounit(theEnergyCutWidth,GeV) << theRapidityCutWidth;
}

void JetRegion::persistentInput(PersistentIStream & is, int) {
  is >> iunit(thePtMin,GeV) >> iunit(thePtMax,GeV)
     >> theYRanges >> theAccepts >> theFuzzy >> theCutWeight
     >> iunit(theEnergyCutWidth,GeV) >> theRapidityCutWidth;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<JetRegion,HandlerBase>
  describeThePEGJetRegion("ThePEG::JetRegion", "JetCuts.so");

void JetRegion::Init() {

  static ClassDocumentation<JetRegion> documentation
    ("JetRegion implements the requirement of finding a jet inside a "
     "given range of transverse momenta, and (pseudo-)rapidity.");

  static Parameter<JetRegion,Energy> interfacePtMin
    ("PtMin",
     "The minimum pt required.",
     &JetRegion::thePtMin, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<JetRegion,Energy> interfacePtMax
    ("PtMax",
     "The maximum pt allowed.",
     &JetRegion::thePtMax, GeV, Constants::MaxEnergy, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Command<JetRegion> interfaceYRange
    ("YRange",
     "Insert a rapidity range.",
     &JetRegion::doYRange, false);

  static ParVector<JetRegion,int> interfaceAccepts
    ("Accepts",
     "The jet numbers accepted. If empty, any jets are accepted.",
     &JetRegion::theAccepts, -1, 1, 1, 10,
     false, false, Interface::upperlim);

  static Switch<JetRegion,bool> interfaceFuzzy
    ("Fuzzy",
     "Make this jet region a fuzzy cut",
     &JetRegion::theFuzzy, false, false, false);
  static SwitchOption interfaceFuzzyOn
    (interfaceFuzzy,
     "Yes",
     "",
     true);
  static SwitchOption interfaceFuzzyOff
    (interfaceFuzzy,
     "No",
     "",
     false);

  static Parameter<JetRegion,Energy> interfaceEnergyCutWidth
    ("EnergyCutWidth",
     "The pt cut smearing.",
     &JetRegion::theEnergyCutWidth, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<JetRegion,double> interfaceRapidityCutWidth
    ("RapidityCutWidth",
     "The rapidity cut smearing.",
     &JetRegion::theRapidityCutWidth, 0.1, 0.0, 0,
     false, false, Interface::lowerlim);

}

