// -*- C++ -*-
//
// DalitzDecayer.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzDecayer class.
//

#include "DalitzDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

DalitzDecayer::~DalitzDecayer() {}

IBPtr DalitzDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr DalitzDecayer::fullclone() const {
  return new_ptr(*this);
}

void DalitzDecayer::doinit() {
  Decayer::doinit();
  rho = getParticleData(ParticleID::rho0);
}

bool DalitzDecayer::accept(const DecayMode & dm) const {
  if ( dm.products().size() != 3 || !dm.cascadeProducts().empty() ||
       !dm.productMatchers().empty() || dm.wildProductMatcher() ) return false;
  bool ep = false, em = false, gam = false;
  for ( ParticleMSet::const_iterator pit = dm.products().begin();
	pit != dm.products().end(); ++pit ) {
    if ( (**pit).id() == ParticleID::eplus ) ep = true;
    else if ( (**pit).id() == ParticleID::eminus ) em = true;
    else if ( (**pit).id() == ParticleID::gamma ) gam = true;
  }
  return ep && em && gam;

}

ParticleVector DalitzDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector children = getChildren(dm, parent);
  tPPtr ep;
  tPPtr em;
  tPPtr gam;
  for ( int i = 0, N = children.size(); i < N; ++i ) {
    if ( children[i]->id() == ParticleID::eplus ) ep = children[i];
    else if ( children[i]->id() == ParticleID::eminus ) em = children[i];
    else if ( children[i]->id() == ParticleID::gamma ) gam = children[i];
  }

  Energy2 mee2 = ZERO;
  Energy2 me2 = ep->mass()*em->mass();
  Energy2 mm2 = sqr(parent.mass());
  Energy2 mr2 = sqr(rho->mass());
  Energy2 gr2 = sqr(rho->width());
  Energy2 mee2min = 4.0*me2;
  do {
    mee2 = mee2min*pow(mm2/mee2min, rnd());
  } while ( rnd() > (1.0 - 2.0*me2/mee2)*sqrt(max(0.0, 1.0 - mee2min/mee2))*
	    pow(1.0 - mee2/mm2, 3.0)*
	    (1.0 + gr2/mr2)/(sqr(1.0 - mee2/mr2) + gr2/mr2) );

  LorentzMomentum pee, p0, pp, pm;
  do {
    SimplePhaseSpace::CMS(mee2, ep, em);
    pee = ep->momentum() + em->momentum();
    p0 = gam->momentum();
    SimplePhaseSpace::CMS(mm2, p0, pee);
    LorentzRotation r (0.0, 0.0, pee.rho()/pee.e());
    r.rotateY(pee.theta());
    r.rotateZ(pee.phi());
    ep->transform(r);
    em->transform(r);
    gam->setMomentum(p0);
    pp = ep->momentum();
    pm = em->momentum();
  } while ( rnd() > ((mee2 - 2.0*me2)*(sqr(p0*pp) + sqr(p0*pm)) +
		     mee2min*((p0*pp)*(p0*pm) + sqr(p0*pp) + sqr(p0*pm)))*4.0/
	    (mee2*sqr(mm2 - mee2)) );
  finalBoost(parent, children);
  setScales(parent, children);

  return children;
}

void DalitzDecayer::persistentOutput(PersistentOStream & os) const {
  os << rho;
}

void DalitzDecayer::persistentInput(PersistentIStream & is, int) {
  is >> rho;
}

void DalitzDecayer::rebind(const TranslationMap & trans)
  {
  rho = trans.translate(rho);
  Decayer::rebind(trans);
}

IVector DalitzDecayer::getReferences() {
  IVector ret = Decayer::getReferences();
  ret.push_back(rho);
  return ret;
}

ClassDescription<DalitzDecayer> DalitzDecayer::initDalitzDecayer;
// Definition of the static class description member.

void DalitzDecayer::Init() {

  static ClassDocumentation<DalitzDecayer> documentation
    ("This class performs Dalitz decays into gamma e+ e-.");

}

