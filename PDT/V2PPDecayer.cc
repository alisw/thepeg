// -*- C++ -*-
//
// V2PPDecayer.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the V2PPDecayer class.
//

#include "V2PPDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

V2PPDecayer::~V2PPDecayer() {}

IBPtr V2PPDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr V2PPDecayer::fullclone() const {
  return new_ptr(*this);
}

bool V2PPDecayer::accept(const DecayMode & dm) const {
  if ( dm.products().size() != 2 || !dm.cascadeProducts().empty() ||
       !dm.productMatchers().empty() || dm.wildProductMatcher() ) return false;
  for ( ParticleMSet::const_iterator it = dm.products().begin();
	it != dm.products().end(); ++it )
    if ( !PseudoScalarMesonMatcher::Check(**it) ) return false;
  return ( VectorMesonMatcher::Check(*dm.parent()) );
}

ParticleVector V2PPDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  grandParent = tPPtr();
  sibling = tPPtr();
  if ( parent.parents().size() == 1 &&
       PseudoScalarMesonMatcher::Check(parent.parents()[0]->data()) )
    grandParent = parent.parents()[0];
  if ( grandParent && grandParent->children().size() == 2 ) {
    tParticleSet siblings = parent.siblings();
    if ( siblings.size() == 1 &&
	 (PseudoScalarMesonMatcher::Check((**siblings.begin()).data()) ||
	  (**siblings.begin()).id() == ParticleID::gamma ) )
      sibling = *siblings.begin();
  }
  return FlatDecayer::decay(dm, parent);
}

double V2PPDecayer::reweight(const DecayMode &, const Particle & parent,
			     const ParticleVector & children) const {
  if ( !sibling || !grandParent ) return 1.0;
  LorentzMomentum gp = grandParent->momentum();
  gp.boost(-parent.momentum().boostVector());
  LorentzMomentum pp(ZERO, ZERO, ZERO, parent.mass());
  Energy2 p10 = pp*gp;
  Energy2 p12 = pp*children[0]->momentum();
  Energy2 p02 = gp*children[0]->momentum();
  Energy2 m02 = gp.m2();
  Energy2 m12 = pp.m2();
  Energy2 m22 = children[0]->momentum().mass2();
  if ( grandParent->id() == ParticleID::gamma )
    return m12*(2.0*p10*p12*p02 - m12*sqr(p02) - m02*sqr(p12) - m22*sqr(p10)
		+ m12*m02*m22)/((sqr(p10) - m12*m02)*(sqr(p12) - m12*m22));
  else
    return sqr(p10*p12 - m12*p02)/((sqr(p10) - m12*m02)*(sqr(p12) - m12*m22));
}

void V2PPDecayer::persistentOutput(PersistentOStream & os) const {
  os << grandParent << sibling;
}

void V2PPDecayer::persistentInput(PersistentIStream & is, int) {
  is >> grandParent >> sibling;
}

ClassDescription<V2PPDecayer> V2PPDecayer::initV2PPDecayer;
// Definition of the static class description member.

void V2PPDecayer::Init() {

  static ClassDocumentation<V2PPDecayer> documentation
    ("This class performs the decay of a vector meson into two "
     "pseudo-scalars according to a flat phase space. If, however the "
     "decaying particle comes from a pseudo-scalar and has only one "
     "sibling which is a pseudo-scalar (or a photon) the decay is "
     "reweighted with cos^2 (sin^2 for photon) of the angle between one "
     "of the decay products and its grand parent. ");

}

