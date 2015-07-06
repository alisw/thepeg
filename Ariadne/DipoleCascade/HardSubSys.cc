// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardSubSys class.
//

#include "HardSubSys.h"
#include "HardRemnant.h"
#include "Ariadne/DipoleCascade/Parton.h"
#include "Ariadne/DipoleCascade/PartonTraits.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HardSubSys.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne;

HardSubSys::~HardSubSys() {}

void HardSubSys::add(tPPtr p) {
  theInitialParticles.push_back(p);
  theInitMomentum = Utilities::sumMomentum(theInitialParticles);
  isModified = true;
}

void HardSubSys::transform(const LorentzRotation & r) {
  isModified = true;
  Utilities::transform(theActivePartons, r);
  Utilities::transform(theProducedParticles, r);
  theInitMomentum.transform(r);
  theRotation.transform(r);
  isRotated = true;
}

void HardSubSys::setMomentum(const LorentzMomentum & q, bool posZDir,
    double x) {
  LorentzMomentum k1, k2, k3, k4;
  //Assume hardRemnants().first is in positive z-direction.
  if(posZDir){
    Energy2 Q2 = hardRemnants().second ? hardRemnants().second->Q2() : 0.0*GeV2;
    Energy pm = momentum().minus() - momentum().perp2() / momentum().plus();
    k1 = lightCone(-Q2 / pm, pm, 0.0*GeV, 0.0*GeV);
  }
  else{
    Energy2 Q2 = hardRemnants().first ? hardRemnants().first->Q2() : 0.0*GeV2;
    Energy pp = momentum().plus() - momentum().perp2() / momentum().minus();
    k1 = lightCone(pp, -Q2/pp, 0.0*GeV, 0.0*GeV);
  }
  //***********************
  k2 = momentum() - k1;
  k3 = k1 + (1 - x) * (q - momentum());
  k4 = k2 + x * (q - momentum());
  //transform(Utilities::getTransformToMomentum(momentum(), q, k));
  transform(Utilities::getBoostFromCM(make_pair(k3, k4))*
      Utilities::getBoostToCM(make_pair(k1, k2)));
  //***********************
  if ( active().size() == 1 && produced().empty() && initial().empty() )
    (**active().begin()).momentum() = q;
}


void HardSubSys::sumMomentum() const {
  if ( !isModified ) return;
  LorentzMomentum tot = theInitMomentum;
  tot += Utilities::sumMomentum(theActivePartons);
  tot += Utilities::sumMomentum(theProducedParticles);
  theMomentum = tot;
  theMomentum.rescaleMass();
  isModified = false;
}

Energy2 HardSubSys::Q2(){
  Energy2 q1  = hardRemnants().first ? hardRemnants().first->Q2() : 0.0*GeV2;
  Energy2 q2  = hardRemnants().second ? hardRemnants().second->Q2() : 0.0*GeV2;
  return q1 + q2;
}

ClonePtr HardSubSys::clone() const {
  return new_ptr(*this);
}

void HardSubSys::fillReferences(CloneSet & cset) const {
  // *** ATTENTION ***
}

void HardSubSys::rebind(const CloneBase::TranslationMap & trans) {
  PartonSet old;
  old.swap(theActivePartons);
  theActivePartons.clear();
  trans.translate(inserter(theActivePartons), old.begin(), old.end());
  old.swap(theProducedParticles);
  theProducedParticles.clear();
  trans.translate(inserter(theProducedParticles), old.begin(), old.end());
  theHardRemnants = make_pair(trans.translate(theHardRemnants.first), 
      trans.translate(theHardRemnants.second));
}

void HardSubSys::persistentOutput(PersistentOStream & os) const {
  os << theActivePartons << theProducedParticles << theIntermediateParticles
     << ounit(theMomentum, GeV) << isModified << ounit(theInitMomentum, GeV)
     << isRotated;
}

void HardSubSys::persistentInput(PersistentIStream & is, int) {
  is >> theActivePartons >> theProducedParticles >> theIntermediateParticles
     >> iunit(theMomentum, GeV) >> isModified >> iunit(theInitMomentum, GeV)
     >> isRotated;
  theRotation =
    Utilities::transformToMomentum(Utilities::sumMomentum(theInitialParticles),
				   theInitMomentum);
}

ClassDescription<HardSubSys> HardSubSys::initHardSubSys;
// Definition of the static class description member.

void HardSubSys::Init() {

}

