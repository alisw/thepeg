// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SoftRemnant class.
//

#include "SoftRemnant.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/EventRecord/Particle.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SoftRemnant.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne;

ClonePtr SoftRemnant::clone() const {
  return new_ptr(*this);
}

void SoftRemnant::parentData(tcPDPtr x) {
  theParentDataPtr = x;
}

double SoftRemnant::partonx() const {
  if( momentum().z() > 0.0*GeV ) {
    return 1.0 - momentum().plus()/parentMomentum().plus();
  }
  else{
    return 1.0 - momentum().minus()/parentMomentum().minus();
  }
}

tPPtr SoftRemnant::produceParticle(const LorentzRotation & r) {
  RemnantParton::produceParticle(LorentzRotation());
  particle()->set5Momentum(r*incomingMomentum());
  return particle();
}

void SoftRemnant::fillReferences(CloneSet & cset) const {
  RemnantParton::fillReferences(cset);
}

void SoftRemnant::rebind(const TranslationMap & trans) {
  RemnantParton::rebind(trans);
}

void SoftRemnant::persistentOutput(PersistentOStream & os) const {
  os << theParentDataPtr << ounit(theParentMomentum, GeV);
}

void SoftRemnant::persistentInput(PersistentIStream & is, int) {
  is >> theParentDataPtr >> iunit(theParentMomentum, GeV);
}

ClassDescription<SoftRemnant> SoftRemnant::initSoftRemnant;
// Definition of the static class description member.

void SoftRemnant::Init() {}

void SoftRemnant::debugme() const {
  RemnantParton::debugme();
  cerr << " sr";
}

