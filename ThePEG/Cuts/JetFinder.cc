// -*- C++ -*-
//
// JetFinder.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2007 Leif Lonnblad
// Copyright (C) 2009-2011 Simon Platzer
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the JetFinder class.
//

#include "JetFinder.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

JetFinder::JetFinder() 
  : theMinOutgoing(1), theRestrictConstituents(false),
    theConstituentRapidityRange(Constants::MaxRapidity,Constants::MaxRapidity) {}

JetFinder::~JetFinder() {}

string JetFinder::doYRange(string in) {
  istringstream ins(in);
  double first, second;
  ins >> first >> second;
  if ( first > second )
    swap(first,second);
  theConstituentRapidityRange = make_pair(first,second);
  return "";
}



// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void JetFinder::persistentOutput(PersistentOStream & os) const {
  os << theUnresolvedMatcher << theMinOutgoing
     << theRestrictConstituents << theConstituentRapidityRange;
}

void JetFinder::persistentInput(PersistentIStream & is, int) {
  is >> theUnresolvedMatcher >> theMinOutgoing 
     >> theRestrictConstituents >> theConstituentRapidityRange;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<JetFinder,Interfaced>
describeJetFinder("ThePEG::JetFinder", "");

void JetFinder::Init() {

  static ClassDocumentation<JetFinder> documentation
    ("JetFinder defines an interface to jet finders to be used when cuts "
     "should actually be defined on the level of reconstructed jets such "
     "as typically encountered in higher order corrections.");


  static Reference<JetFinder,MatcherBase> interfaceUnresolvedMatcher
    ("UnresolvedMatcher",
     "A matcher identifying unresolved partons",
     &JetFinder::theUnresolvedMatcher, false, false, true, false, false);


  static Parameter<JetFinder,unsigned int> interfaceMinOutgoing
    ("MinOutgoing",
     "The minimum number of outgoing partons to be clustered.",
     &JetFinder::theMinOutgoing, 1, 1, 0,
     false, false, Interface::lowerlim);


  static Switch<JetFinder,bool> interfaceRestrictConstituents
    ("RestrictConstituents",
     "Restrict the constituents for clustering.",
     &JetFinder::theRestrictConstituents, false, false, false);
  static SwitchOption interfaceRestrictConstituentsYes
    (interfaceRestrictConstituents,
     "Yes",
     "",
     true);
  static SwitchOption interfaceRestrictConstituentsNo
    (interfaceRestrictConstituents,
     "No",
     "",
     false);

  static Command<JetFinder> interfaceYRange
    ("ConstituentRapidityRange",
     "Restrict clustering to a rapidity range.",
     &JetFinder::doYRange, false);

}

