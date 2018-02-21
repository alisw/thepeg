// -*- C++ -*-
//
// Tau2HadronsDecayer.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Tau2HadronsDecayer class.
//

#include "Tau2HadronsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

Tau2HadronsDecayer::~Tau2HadronsDecayer() {}

IBPtr Tau2HadronsDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr Tau2HadronsDecayer::fullclone() const {
  return new_ptr(*this);
}

bool Tau2HadronsDecayer::accept(const DecayMode & dm) const {
  if ( dm.products().size() < 3 || !dm.cascadeProducts().empty() ||
       !dm.productMatchers().empty() || dm.wildProductMatcher() ) return false;
  if ( abs(dm.parent()->id()) != ParticleID::tauminus ) return false;
  for ( ParticleMSet::const_iterator pit = dm.products().begin();
	pit != dm.products().end(); ++pit )
    if ( (**pit).id()*dm.parent()->id() ==
	 ParticleID::tauminus*ParticleID::nu_tau ) return true;
  return false;
}

double Tau2HadronsDecayer::reweight(const DecayMode &, const Particle & parent,
				    const ParticleVector & children) const {
  tPPtr nu;
  for ( int i = 0, N = children.size(); i < N; ++i )
    if ( parent.id()*children[i]->id() ==
	 ParticleID::tauminus*ParticleID::nu_tau ) {
      nu = children[i];
      break;
    }
  double xnu = 2.0*nu->momentum().e()/parent.mass();
  return 0.5*xnu*(3.0 - xnu);
}

void Tau2HadronsDecayer::persistentOutput(PersistentOStream &) const {}

void Tau2HadronsDecayer::persistentInput(PersistentIStream &, int) {}

ClassDescription<Tau2HadronsDecayer> Tau2HadronsDecayer::initTau2HadronsDecayer;
// Definition of the static class description member.

void Tau2HadronsDecayer::Init() {

  static ClassDocumentation<Tau2HadronsDecayer> documentation
    ("This class will perform the decays of tau to neutrimo + hadrons "
      "according to phase space, with an extra weight xnu(3-nxu).");

}

