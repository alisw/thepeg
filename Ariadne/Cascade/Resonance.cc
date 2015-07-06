// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Resonance class.
//

#include "Resonance.h"
#include "Emission.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

Resonance::Resonance() {}

Resonance::~Resonance() {}

ClonePtr Resonance::clone() const {
  return new_ptr(*this);
}

void Resonance::fillReferences(CloneSet & cset) const {
  Parton::fillReferences(cset);
}

void Resonance::rebind(const TranslationMap & trans) {
  Parton::rebind(trans);
  theParentResonance = trans.translate(theParentResonance);
}

void Resonance::persistentOutput(PersistentOStream & os) const {
  os << theParentResonance << theDecaySystem;
}

void Resonance::persistentInput(PersistentIStream & is, int) {
  is >> theParentResonance >> theDecaySystem;
}


// The following static variable is needed for the type description
// system in ThePEG.
DescribeClass<Resonance,Parton>
  describeAriadne5Resonance("Ariadne5::Resonance", "libAriadne5.so");

void Resonance::Init() {}

