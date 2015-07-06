// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ResonanceFinder class.
//

#include "ResonanceFinder.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Config/algorithm.h"



using namespace Ariadne5;

ResonanceFinder::ResonanceFinder() {}

ResonanceFinder::~ResonanceFinder() {}

IBPtr ResonanceFinder::clone() const {
  return new_ptr(*this);
}

IBPtr ResonanceFinder::fullclone() const {
  return new_ptr(*this);
}


tPVector ResonanceFinder::resonances(SubProcess & sub) const {
  list<tPPtr> found;
  for ( int i = 0, N = sub.intermediates().size(); i < N; ++i ) {
    tPPtr res = sub.intermediates()[i];
    if ( res->mass() <= ZERO ) continue;
    list<tPPtr>::iterator it = found.begin();
    while ( it != found.end() && !member(res->children(), *it) ) ++it;
    found.insert(it, res);
  }
  return tPVector(found.begin(), found.end());
}

void ResonanceFinder::insert(Step & step, SubProcess & sub, PPtr resonance,
			     const tPVector & children) const {
  for ( int i = 0, N = children.size(); i < N; ++i ) {
    for ( int j = 0, M = children[i]->parents().size(); j < M; ++j ) {
      if ( !member(children, children[i]->parents()[j]) ) {
	children[i]->parents()[j]->abandonChild(children[i]);
	children[i]->parents()[j]->addChild(resonance);
      }
    }
    resonance->addChild(children[i]);
  }
  sub.addIntermediate(resonance, false);
  step.addIntermediate(resonance);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

// The following static variable is needed for the type description
// system in ThePEG.
DescribeNoPIOClass<ResonanceFinder,Interfaced>
  describeAriadne5ResonanceFinder("Ariadne5::ResonanceFinder", "libAriadne5.so");

void ResonanceFinder::Init() {

  static ClassDocumentation<ResonanceFinder> documentation
    ("The ResonanceFinder class is used by the AriadneHandler to find "
     "possible s-channel resonances in a given SubProcess. Normally this "
     "is quite trivial, but in some cases the information about resonances "
     "is not available in the SubProcess (eg. if the process has been "
     "read in from an Les Houches event file). In this case the user may "
     "inherit from this class to specify which of the final state "
     "particles in the SubProcess was produced via a resonance.");

}

