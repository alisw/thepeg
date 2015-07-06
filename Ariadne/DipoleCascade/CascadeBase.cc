// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CascadeBase class.
//

#include "CascadeBase.h"
#include "DipoleState.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CascadeBase.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne;

void CascadeBase::fillReferences(CloneSet & cset) const {
  cset.insert(theState);
}

void CascadeBase::rebind(const TranslationMap & trans) {
  theState = trans.translate(theState);
}

void CascadeBase::persistentOutput(PersistentOStream & os) const {
  os << theHandler << theState << isTouched;
}

void CascadeBase::persistentInput(PersistentIStream & is, int) {
  is >> theHandler >> theState >> isTouched;
}

AbstractClassDescription<CascadeBase> CascadeBase::initCascadeBase;
// Definition of the static class description member.

void CascadeBase::Init() {}

