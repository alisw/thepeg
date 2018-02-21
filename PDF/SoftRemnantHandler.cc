// -*- C++ -*-
//
// SoftRemnantHandler.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SoftRemnantHandler class.
//

#include "SoftRemnantHandler.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Utilities/Direction.h"
#include "ThePEG/Utilities/Maths.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/RemnantDecayer.h"
#include "ThePEG/PDT/RemnantData.h"
#include "ThePEG/EventRecord/RemnantParticle.h"

using namespace ThePEG;

IBPtr SoftRemnantHandler::clone() const {
  return new_ptr(*this);
}

IBPtr SoftRemnantHandler::fullclone() const {
  return new_ptr(*this);
}

bool SoftRemnantHandler::
canHandle(tcPDPtr particle, const cPDVector & partons) const {
  if ( !remdec ) return false;
  for ( int i = 0, N = partons.size(); i < N; ++i )
    if ( !remdec->canHandle(particle, partons[i]) ) return false;
  return true;
}

Lorentz5Momentum SoftRemnantHandler::
generate(PartonBinInstance & pb, const double *,
	 Energy2, const LorentzMomentum & parent,
	 bool fixedPartonMomentum) const {
  if ( !fixedPartonMomentum ) {
    LorentzMomentum p = lightCone((parent.rho() + parent.e())*pb.xi(), Energy());
    p.rotateY(parent.theta());
    p.rotateZ(parent.phi());
    pb.parton()->setMomentum(p);
  }
  RemPPtr rem = new_ptr(RemnantParticle(*pb.particle(), remdec, pb.parton()));
  if ( rem->extracted().empty() ) pb.remnantWeight(0.0);
  pb.remnants(PVector(1, rem));
  return pb.parton()->momentum();
}

bool SoftRemnantHandler::
recreateRemnants(PartonBinInstance & pb, tPPtr oldp, tPPtr newp, double,
		 Energy2, const LorentzMomentum &,
		 const PVector &) const {
  // First find the old remnant.
  RemPPtr rem;
  for ( int i = 0, N = pb.particle()->children().size(); i < N; ++i )
    if ( dynamic_ptr_cast<tRemPPtr>(pb.particle()->children()[i]) )
      rem = dynamic_ptr_cast<tRemPPtr>(pb.particle()->children()[i]);
  if ( !rem ) return false;

//   LorentzMomentum p(0.0, 0.0, parent.rho(), parent.e());
//   pb.parton()->setMomentum
//     (lightCone(p.plus()*pb.x(), ZERO, ZERO, ZERO));
//   pb.parton()->rotateY(parent.theta());
//   pb.parton()->rotateZ(parent.phi());

  rem = new_ptr(*rem);
  if ( !rem->reextract(oldp, newp) ) return false;

  pb.remnants(PVector(1, rem));

  return true;
}  

Lorentz5Momentum SoftRemnantHandler::
generate(PartonBinInstance & pb, const double *, Energy2, Energy2,
	 const LorentzMomentum & parent, bool fixedPartonMomentum) const {
  return generate(pb, 0, ZERO, parent, fixedPartonMomentum);
}

bool SoftRemnantHandler::
recreateRemnants(PartonBinInstance & pb, tPPtr oldp, tPPtr newp, double l,
		 Energy2 scale, Energy2,
		 const LorentzMomentum & p, const PVector & prev) const {
  return recreateRemnants(pb, oldp, newp, l, scale, p, prev);
}  

void SoftRemnantHandler::persistentOutput(PersistentOStream & os) const {
  os << remdec;
}

void SoftRemnantHandler::persistentInput(PersistentIStream & is, int) {
  is >> remdec;
}

void SoftRemnantHandler::setDecayer(RemDecPtr rd) {
  remdec = rd;
  isMultiCapable = rd->multiCapable();
}

ClassDescription<SoftRemnantHandler>
SoftRemnantHandler::initSoftRemnantHandler;

void SoftRemnantHandler::Init() {

  static ClassDocumentation<SoftRemnantHandler> documentation
    ("SoftRemnantHandler inherits from the RemnantHandler and implements "
     "the generation of a single collinear RemnantParticle when anything "
     "is extracted from anything else. Such a RemnantParticle needs to be "
     "decayed by a special RemnantDecayer and the SoftRemnantHandler "
     "needs to be assign such a decayer to work properly.");

  static Reference<SoftRemnantHandler,RemnantDecayer> interfaceRemnantDecayer
    ("RemnantDecayer",
     "A RemnantDecayer object which is able to decay the produced "
     "RemnantParticle objects.",
     &SoftRemnantHandler::remdec, false, false, true, false, true,
     &SoftRemnantHandler::setDecayer);

}

