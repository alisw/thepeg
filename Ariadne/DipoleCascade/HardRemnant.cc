// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardRemnant class.
//

#include "HardRemnant.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HardRemnant.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne;

using namespace ThePEG;

HardRemnant::~HardRemnant() {}

ClonePtr HardRemnant::clone() const {
  return new_ptr(*this);
}

void HardRemnant::fillReferences(CloneSet & cset) const {
  RemnantParton::fillReferences(cset);
}

void HardRemnant::rebind(const TranslationMap & trans) {
  RemnantParton::rebind(trans);
}

void HardRemnant::persistentOutput(PersistentOStream & os) const {
  os << ounit(theQ2, GeV2);
}

void HardRemnant::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theQ2, GeV2);
}

ClassDescription<HardRemnant> HardRemnant::initHardRemnant;
// Definition of the static class description member.

void HardRemnant::Init() {}

void HardRemnant::debugme() const {
  RemnantParton::debugme();
  cerr << " hr";
}

