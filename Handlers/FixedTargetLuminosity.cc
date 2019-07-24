// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FixedTargetLuminosity class.
//

#include "FixedTargetLuminosity.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

FixedTargetLuminosity::FixedTargetLuminosity() {}

FixedTargetLuminosity::~FixedTargetLuminosity() {}

IBPtr FixedTargetLuminosity::clone() const {
  return new_ptr(*this);
}

IBPtr FixedTargetLuminosity::fullclone() const {
  return new_ptr(*this);
}

void FixedTargetLuminosity::persistentOutput(PersistentOStream & os) const {
  os << beam_ << target_ << ounit(ecms_,GeV) << beta_;
}

void FixedTargetLuminosity::persistentInput(PersistentIStream & is, int) {
  is >> beam_ >> target_ >> iunit(ecms_,GeV) >> beta_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FixedTargetLuminosity,LuminosityFunction>
describeThePEGFixedTargetLuminosity("ThePEG::FixedTargetLuminosity",
				    "FixedTargetLuminosity.so");

void FixedTargetLuminosity::Init() {

  static ClassDocumentation<FixedTargetLuminosity> documentation
    ("The FixedTargetLuminosity class implements the luminosity for"
     "fixed target collisions");

  static Reference<FixedTargetLuminosity,ParticleData> interfaceTargetParticle
    ("TargetParticle",
     "The target particle",
     &FixedTargetLuminosity::target_, false, false, true, false, false);

  static Reference<FixedTargetLuminosity,ParticleData> interfaceBeamParticle
    ("BeamParticle",
     "The beam particle",
     &FixedTargetLuminosity::beam_, false, false, true, false, false);

}

bool FixedTargetLuminosity::canHandle(const cPDPair & pdpair) const {
  return pdpair.first==beam_ && pdpair.second==target_;
}

Energy FixedTargetLuminosity::maximumCMEnergy() const {
  return ecms_;
}

LorentzRotation FixedTargetLuminosity::getBoost() const {
  return LorentzRotation(0.0, 0.0, beta_);
}

double FixedTargetLuminosity::Y() const {
  Energy en = beamEMaxA()+target_->mass();
  Energy pp = sqrt(sqr(beamEMaxA())+sqr(beam_->mass())); 
  return 0.5*log((en+pp)/(en-pp));
}

void FixedTargetLuminosity::doinit() {
  Energy2 m12 = sqr(beam_  ->mass());
  Energy2 m22 = sqr(target_->mass());
  Energy2 Em2 = beamEMaxA()*target_->mass();
  Energy2 p2 = (sqr(Em2)-m12*m22)/(2.*Em2+m12+m22);
  beta_ = sqrt(p2/(p2+m22));
  ecms_ = sqrt(2.*Em2+m12+m22);
  LuminosityFunction::doinit();
}
