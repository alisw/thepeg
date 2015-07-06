// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RemnantParton class.
//

#include "RemnantParton.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/EventRecord/Particle.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "RemnantParton.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne;

ClonePtr RemnantParton::clone() const {
  return new_ptr(*this);
}

void RemnantParton::fillReferences(CloneSet & cset) const {
  Parton::fillReferences(cset);
}

void RemnantParton::rebind(const TranslationMap & trans) {
  Parton::rebind(trans);
}

void RemnantParton::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMu, GeV) << theAlpha;
}

void RemnantParton::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMu, GeV) >> theAlpha;
}

ClassDescription<RemnantParton> RemnantParton::initRemnantParton;
// Definition of the static class description member.

void RemnantParton::Init() {}

void RemnantParton::debugme() const {
  Parton::debugme();
}

