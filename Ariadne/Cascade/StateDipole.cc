// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StateDipole class.
//

#include "StateDipole.h"
#include "AriadneHandler.h"

#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Ariadne5;

StateDipole::StateDipole() {
  theRhoCut = Current<AriadneHandler>()->pTCut();
}

ClonePtr StateDipole::clone() const {
  return new_ptr(*this);
}

void StateDipole::fillReferences(CloneSet & cset) const {
  DipoleBase::fillReferences(cset);
}

void StateDipole::rebind(const TranslationMap & trans) {
  DipoleBase::rebind(trans);
}

DescribeNoPIOClass<StateDipole,DipoleBase>
describeAriadne5StateDipole("Ariadne5::StateDipole", "libAriadne5.so");

void StateDipole::Init() {}

