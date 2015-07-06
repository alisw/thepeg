// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Parton class.
//

#include "Parton.h"
#include "DipoleState.h"
#include "Dipole.h"
#include "String.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Particle.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Parton.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne;

Parton::Parton()
  : isGluon(false) {}

Parton::Parton(const Parton & x)
  : CascadeBase(x), theOrig(x.theOrig), theParticle(x.theParticle),
    theParents(x.theParents),
    theDataPtr(x.theDataPtr), isGluon(x.isGluon),
    theMomentum(x.theMomentum), theString(x.theString),
    theIDip(x.theIDip), theODip(x.theODip) {}

Parton::~Parton() {}

ClonePtr Parton::clone() const {
  return new_ptr(*this);
}

void Parton::orig(tcPPtr x) {
  theOrig = x;
  data(x->dataPtr());
  theMomentum = x->momentum();
}

tPPtr Parton::produceParticle(const LorentzRotation & r) {
  theParticle = data().produceParticle(r*momentum());
  particle()->scale(sqr(coloured()? handler()->pTCut(): handler()->pTCutEM()));
  return particle();
}

void Parton::data(tcPDPtr x) {
  theDataPtr = x;
  isGluon = ( data().id() == ParticleID::g );
}

tParPtr Parton::prev() const {
  return iDip()? iDip()->iPart(): tParPtr();
}

tParPtr Parton::next() const {
  return oDip()? oDip()->oPart(): tParPtr();
}

Energy2 Parton::invPT2() const {
  if(!next() || !prev()){
    return 0.0*GeV2;
  }
  Energy2 s123 = (next()->momentum() + prev()->momentum() + momentum()).m2();
  Energy2 s12 = (prev()->momentum() + momentum()).m2();
  Energy2 s23 = (next()->momentum() + momentum()).m2();
  Energy m1 = prev()->momentum().mass();
  Energy m2 = momentum().mass();
  Energy m3 = next()->momentum().mass();
  return (s12 - sqr(m1 + m2)) * (s23 - sqr(m2 + m3)) / s123;
}

void Parton::fillReferences(CloneSet & cset) const {
  CascadeBase::fillReferences(cset);
  cset.insert(iDip());
  cset.insert(oDip());
}

void Parton::rebind(const TranslationMap & trans) {
  CascadeBase::rebind(trans);
  theIDip = trans.translate(theIDip);
  theODip = trans.translate(theODip);
  theParents.first = trans.translate(theParents.first);
  theParents.second = trans.translate(theParents.second);
  theString = trans.translate(theString);
}

void Parton::persistentOutput(PersistentOStream & os) const {
  os << theOrig << theParticle << theParents << theDataPtr << isGluon
     << ounit(theMomentum, GeV) << theString << theIDip << theODip;
}

void Parton::persistentInput(PersistentIStream & is, int) {
  is >> theOrig >> theParticle >> theParents >> theDataPtr >> isGluon
     >> iunit(theMomentum, GeV) >> theString >> theIDip >> theODip;
}

ClassDescription<Parton> Parton::initParton;
// Definition of the static class description member.

void Parton::Init() {}

void Parton::debugme() const {
  CascadeBase::debugme();
  cerr << "P"
       << setw(3) << state()->index(this)
       << setw(4) << data().id()
       << setw(3) << state()->index(string())
       << setw(3) << state()->index(iDip())
       << setw(3) << state()->index(oDip())
       << (touched()? " *": "  ") << setprecision(3)
       << setw(9)
       << ( abs(momentum().x()/momentum().e()) < 1e-6?
	    0.0*GeV: momentum().x())/GeV
       << setw(9)
       << ( abs(momentum().y()/momentum().e()) < 1e-6?
	    0.0*GeV: momentum().y())/GeV
       << setw(9)
       << ( abs(momentum().z()/momentum().e()) < 1e-6?
	    0.0*GeV: momentum().z())/GeV
       << setw(9) << momentum().e()/GeV
       << setw(9)
       << ( abs(momentum().m()/momentum().e()) < 1e-6?
	    0.0*GeV: momentum().m())/GeV
       << setw(9) << momentum().mass()/GeV;
}

