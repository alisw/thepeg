// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CascadeBase class.
//

#include "CascadeBase.h"
#include "DipoleState.h"

#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

void CascadeBase::fillReferences(CloneSet & cset) const {}

void CascadeBase::rebind(const TranslationMap & trans) {
  theState = trans.translate(theState);
}

void CascadeBase::persistentOutput(PersistentOStream & os) const {
  os << theState << isTouched;
}

void CascadeBase::persistentInput(PersistentIStream & is, int) {
  is >> theState >> isTouched;
}

// Definition of the static class description member.
DescribeAbstractClass<CascadeBase,CloneBase>
describeAriadne5CascadBase("Ariadne5::CascadeBase", "libAriadne5.so");

void CascadeBase::Init() {}

