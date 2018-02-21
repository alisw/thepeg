// -*- C++ -*-
//
// OmegaPhi3PiDecayer.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OmegaPhi3PiDecayer class.
//

#include "OmegaPhi3PiDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

IBPtr OmegaPhi3PiDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr OmegaPhi3PiDecayer::fullclone() const {
  return new_ptr(*this);
}

bool OmegaPhi3PiDecayer::accept(const DecayMode & dm) const {
  if ( dm.products().size() != 3 || !dm.cascadeProducts().empty() ||
       !dm.productMatchers().empty() || dm.wildProductMatcher() ) return false;
  bool foundp = false, foundm = false, found0 = false;
  for ( ParticleMSet::const_iterator pit = dm.products().begin();
	pit != dm.products().end(); ++pit ) {
    if ( (**pit).id() == ParticleID::piplus ) foundp = true;
    else if ( (**pit).id() == ParticleID::piminus ) foundm = true;
    else if ( (**pit).id() == ParticleID::pi0 ) found0 = true;
  }
  return foundp && foundm && found0;
}

double OmegaPhi3PiDecayer::reweight(const DecayMode &, const Particle & parent,
				    const ParticleVector & children) const {
  Energy2 m2 = parent.momentum().mass2();
  Lorentz5Momentum pp;
  Lorentz5Momentum pm;
  Lorentz5Momentum p0;
  for ( int i = 0, N = children.size(); i < N; ++i ) {
    if ( children[i]->id() == ParticleID::piplus )
      pp = children[i]->momentum();
    else if ( children[i]->id() == ParticleID::piminus )
      pm = children[i]->momentum();
    else if ( children[i]->id() == ParticleID::pi0 )
      p0 = children[i]->momentum();
  }

  return (pp.mass2()*pm.mass2()*p0.mass2() + 2.0*(pp*pm)*(pp*p0)*(pm*p0)
	  - pp.mass2()*sqr(pm*p0) - pm.mass2()*sqr(pp*p0)
	  - p0.mass2()*sqr(pp*pm))*margin/(m2*m2*m2);
}

void OmegaPhi3PiDecayer::persistentOutput(PersistentOStream & os) const {
  os << margin;
}

void OmegaPhi3PiDecayer::persistentInput(PersistentIStream & is, int) {
  is >> margin;
}

ClassDescription<OmegaPhi3PiDecayer> OmegaPhi3PiDecayer::initOmegaPhi3PiDecayer;
// Definition of the static class description member.

void OmegaPhi3PiDecayer::Init() {

  static ClassDocumentation<OmegaPhi3PiDecayer> documentation
    ("This class performs the decay of a phi or an omega into pi+ pi- "
     "p0. It will in fact decay anything into pi+ pi- p0 assuming the "
     "same matrix element.");


  static Parameter<OmegaPhi3PiDecayer,double> interfaceMargin
    ("Margin",
     "Used to multiply the bare weight to get something below unity. "
     "In the Fortran pythia version it was set to 150 for unknown reasons.",
     &OmegaPhi3PiDecayer::margin, 150.0, 0.0, 1000.0,
     true, false, true);


}

