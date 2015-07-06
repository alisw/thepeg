// -*- C++ -*-
//
// UnResolvedRemnant.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UnResolvedRemnant class.
//

#include "UnResolvedRemnant.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Utilities/Direction.h"
#include "ThePEG/Utilities/Maths.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

IBPtr UnResolvedRemnant::clone() const {
  return new_ptr(*this);
}

IBPtr UnResolvedRemnant::fullclone() const {
  return new_ptr(*this);
}

UnResolvedRemnant::UnResolvedRemnant()
  : minX(1.0e-10) {}

void UnResolvedRemnant::doinit() {
  thePhoton = getParticleData(ParticleID::gamma);
  RemnantHandler::doinit();
}

bool UnResolvedRemnant::
canHandle(tcPDPtr particle, const cPDVector & partons) const {
  for ( cPDVector::const_iterator it = partons.begin();
	it != partons.end(); ++it ) {
    if ( (**it).id() != particle->id()    &&
	 ( (**it).id() != ParticleID::gamma &&
	   (**it).id() != ParticleID::pomeron &&
	   (**it).id() != ParticleID::reggeon) ) return false;
  }
  return true;
}

int UnResolvedRemnant::nDim(const PartonBin &, bool) const {
  return 1;
}

Lorentz5Momentum UnResolvedRemnant::
generate(PartonBinInstance & pb, const double *r,
	 Energy2 scale, const LorentzMomentum & parent,
	 bool fixedPartonMomentum) const {
  // photon into hard process and lepton remnant
  if ( pb.particleData() != pb.partonData()) {
    scale = abs(scale);
    Energy  ppl = pb.xi()*(abs(parent.z())+parent.t());
    Energy2 qt2 = pb.eps()*scale-sqr(pb.xi()*parent.m());
    if(qt2<ZERO) {
      pb.remnantWeight(-1.0);
      return Lorentz5Momentum();
    }
    Energy  pmi = (qt2-scale)/ppl;
    Lorentz5Momentum pgam;
    if ( !fixedPartonMomentum ) {
      pgam.setMass(-sqrt(scale));
      pgam.setT(0.5*(ppl+pmi));
      pgam.setZ(0.5*(ppl-pmi));
      double phi = r[0]*Constants::twopi;
      pgam.setX(sqrt(qt2)*cos(phi));
      pgam.setY(sqrt(qt2)*sin(phi));
      pgam.rotateY(parent.theta());
      pgam.rotateZ(parent.phi());
    } else {
      pgam = pb.parton()->momentum();
    }
    Lorentz5Momentum prem=parent-pgam;
    PPtr rem = pb.particleData()->produceParticle(prem, pb.particleData()->mass());
    pb.remnants(PVector(1, rem));
    return pgam;
  }
  else {
    if ( pb.eps() < minX ) {
      pb.remnants(PVector());
      return parent;
    }
    LorentzMomentum p(ZERO, ZERO, parent.rho(), parent.e());
    TransverseMomentum qt;
    Energy2 qt2 = ZERO;
    if ( scale >= ZERO ) {
      qt2 = pb.eps()*(pb.xi()*parent.m2() + scale);
      double phi = r[0]*Constants::twopi;
      qt = TransverseMomentum(sqrt(qt2)*cos(phi), sqrt(qt2)*sin(phi));
    }
    Energy pl = p.plus()*pb.eps();
    LorentzMomentum prem;
    if ( !fixedPartonMomentum ) {
      prem = lightCone(pl, qt2/pl, qt);
      prem.rotateY(parent.theta());
      prem.rotateZ(parent.phi());
    } else {
      prem = parent - pb.parton()->momentum();
    }
    PPtr rem = thePhoton->produceParticle(prem, ZERO);
    pb.remnants(PVector(1, rem));
    return parent - rem->momentum();
  }
}

Lorentz5Momentum UnResolvedRemnant::
generate(PartonBinInstance & pb, const double *r, Energy2 scale, Energy2,
	 const LorentzMomentum & parent,
	 bool fixedPartonMomentum) const {
  // photon into hard process and lepton remnant
  if ( pb.particleData() != pb.partonData()) {
    scale = abs(scale);
    Energy  ppl = pb.xi()*(abs(parent.z())+parent.t());
    Energy2 qt2 = pb.eps()*scale-sqr(pb.xi()*parent.m());
    if(qt2<ZERO) {
      pb.remnantWeight(-1.0);
      return Lorentz5Momentum();
    }
    Energy  pmi = (qt2-scale)/ppl;
    Lorentz5Momentum pgam;
    if ( !fixedPartonMomentum ) {
      pgam.setMass(-sqrt(scale));
      pgam.setT(0.5*(ppl+pmi));
      pgam.setZ(0.5*(ppl-pmi));
      double phi = r[0]*Constants::twopi;
      pgam.setX(sqrt(qt2)*cos(phi));
      pgam.setY(sqrt(qt2)*sin(phi));
      pgam.rotateY(parent.theta());
      pgam.rotateZ(parent.phi());
    } else {
      pgam = pb.parton()->momentum();
    }
    Lorentz5Momentum prem=parent-pgam;
    PPtr rem = pb.particleData()->produceParticle(prem, pb.particleData()->mass());
    pb.remnants(PVector(1, rem));
    return pgam;
  }
  else {
    if ( pb.eps() < minX ) {
      pb.remnants(PVector());
      return parent;
    }
    LorentzMomentum p(ZERO, ZERO, parent.rho(), parent.e());
    TransverseMomentum qt;
    Energy2 qt2 = ZERO;
    if ( scale >= ZERO ) {
      qt2 = pb.eps()*(pb.xi()*parent.m2() + scale);
      double phi = r[0]*Constants::twopi;
      qt = TransverseMomentum(sqrt(qt2)*cos(phi), sqrt(qt2)*sin(phi));
    }
    Energy pl = p.plus()*pb.eps();
    LorentzMomentum prem;
    if ( !fixedPartonMomentum ) {
      prem = lightCone(pl, qt2/pl, qt);
      prem.rotateY(parent.theta());
      prem.rotateZ(parent.phi());
    } else {
      prem = parent - pb.parton()->momentum();
    }
    PPtr rem = thePhoton->produceParticle(prem, ZERO);
    pb.remnants(PVector(1, rem));
    return parent - rem->momentum();
  }
}

void UnResolvedRemnant::persistentOutput(PersistentOStream & os) const {
  os << minX << thePhoton;
}

void UnResolvedRemnant::persistentInput(PersistentIStream & is, int) {
  is >> minX >> thePhoton;
}

ClassDescription<UnResolvedRemnant>
UnResolvedRemnant::initUnResolvedRemnant;

void UnResolvedRemnant::Init() {

  static ClassDocumentation<UnResolvedRemnant> documentation
    ("UnResolvedRemnant inherits from the RemnantHandler and implements"
     "the generation of either the incoming particle as the remnant"
     "with the emission of a photon, pomeron or reggeon, or"
     "a photon remnant for the particle entering the hard process.");

  static Parameter<UnResolvedRemnant,double> interfaceMinX
    ("MinX",
     "The minimum energy fraction allowed for a photon remnant. "
     "If less than this no remnant will be emitted.",
     &UnResolvedRemnant::minX, 1.0e-10, 0.0, 1.0,
     true, false, true);

  interfaceMinX.rank(10);

}


bool UnResolvedRemnant::
recreateRemnants(PartonBinInstance & pb, tPPtr oldp, tPPtr newp, double,
		 Energy2 scale, const LorentzMomentum & p,
		 const PVector & prev) const {
  if ( !oldp || !prev.empty() ) return false;
  // get the same random number used for the azimuth last time
  Lorentz5Momentum pgam=oldp->momentum();
  pgam.rotateZ(-p.phi());
  pgam.rotateY(-p.theta());
  double test = atan2(pgam.y(),pgam.x())/Constants::twopi;
  if(test<0.) test+=1.;
  vector<double> rv;
  int rd = pb.bin()->remDim();
  for ( int i = 0; i < rd; ++i) rv.push_back(test);
  // compute the momentum
  newp->set5Momentum(generate(pb, pb.bin()->remDim()? &rv[0]: 0, scale, p));
  boostRemnants(pb);
  return true;
}  

bool UnResolvedRemnant::
recreateRemnants(PartonBinInstance & pb, tPPtr oldp, tPPtr newp, double,
		 Energy2 scale, Energy2 shat,
		 const LorentzMomentum & p, const PVector & prev) const {
  if ( !oldp || !prev.empty() ) return false;
  // get the same random number used for the azimuth last time
  Lorentz5Momentum pgam=oldp->momentum();
  pgam.rotateZ(-p.phi());
  pgam.rotateY(-p.theta());
  double test = atan2(pgam.y(),pgam.x())/Constants::twopi;
  if(test<0.) test+=1.;
  vector<double> rv;
  int rd = pb.bin()->remDim();
  for ( int i = 0; i < rd; ++i) rv.push_back(test);
  // compute the momentum
  newp->set5Momentum(generate(pb, rd > 0? &rv[0]: 0, scale, shat, p));
  boostRemnants(pb);
  return true;
}  
