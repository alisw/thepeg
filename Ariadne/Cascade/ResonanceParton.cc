// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ResonanceParton class.
//

#include "ResonanceParton.h"
#include "Emission.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

ResonanceParton::ResonanceParton(): Parton(true) {}

ResonanceParton::~ResonanceParton() {}

ClonePtr ResonanceParton::clone() const {
  return new_ptr(*this);
}

void ResonanceParton::fillReferences(CloneSet & cset) const {
  Parton::fillReferences(cset);
}

void ResonanceParton::rebind(const TranslationMap & trans) {
  Parton::rebind(trans);
  theResonance = trans.translate(theResonance);
  set<tParPtr> oldsib;
  oldsib.swap(theSiblings);
  theSiblings.clear();
  trans.translate(inserter(theSiblings), oldsib.begin(), oldsib.end());
  set<tParPtr> oldoth;
  oldoth.swap(theOthers);
  theOthers.clear();
  trans.translate(inserter(theOthers), oldoth.begin(), oldoth.end());
}

void ResonanceParton::notify(const Emission & em) {
  for ( int i = 0, N = em.partons.size(); i < N; ++i )
    if ( em.partons[i]->parents().first &&
	 member(siblings(), em.partons[i]->parents().first) &&
	 em.partons[i]->parents().first &&
	 member(siblings(), em.partons[i]->parents().first) )
      addSibling(em.partons[i]);
  for ( int i = 0, N = em.radiators.size(); i < N; ++i )
    if ( !em.radiators[i]->finalState() &&
	 member(siblings(), em.radiators[i]) )
      addOther(em.radiators[i]);
  for ( int i = 0, N = em.affected.size(); i < N; ++i )
    if ( !em.affected[i]->finalState() && member(siblings(), em.affected[i]) )
      addOther(em.affected[i]);
}

void ResonanceParton::addOther(tParPtr o) {
    theOthers.insert(o);
    theSiblings.erase(o);
}

void ResonanceParton::persistentOutput(PersistentOStream & os) const {
  os << theResonance << theSiblings << theOthers;
}

void ResonanceParton::persistentInput(PersistentIStream & is, int) {
  is >> theResonance >> theSiblings >> theOthers;
}


// The following static variable is needed for the type description
// system in ThePEG.
DescribeClass<ResonanceParton,Parton>
  describeAriadne5ResonanceParton("Ariadne5::ResonanceParton",
				"libAriadne5.so");

void ResonanceParton::Init() {}

