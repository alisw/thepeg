// -*- C++ -*-
//
// Onium3GDecayer.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Onium3GDecayer class.
//

#include "Onium3GDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/HandlerGroup.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

Onium3GDecayer::~Onium3GDecayer() {}

IBPtr Onium3GDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr Onium3GDecayer::fullclone() const {
  return new_ptr(*this);
}

bool Onium3GDecayer::accept(const DecayMode & dm) const {
  if ( dm.products().size() != 3 || !dm.cascadeProducts().empty() ||
       !dm.productMatchers().empty() || dm.wildProductMatcher() ) return false;
  if ( !VectorMesonMatcher::Check(*dm.parent()) ) return false;
  vector<long> flv = PDT::flavourContent(*dm.parent());
  if ( abs(flv[0]) < 4 || flv[0] + flv[1] != 0 ) return false;
  int ng = 0;
  int np = 0;
  for ( int i = 0; i < 3; ++i )
    if ( dm.orderedProducts()[i]->id() == ParticleID::g ) ++ng;
    else if ( dm.orderedProducts()[i]->id() == ParticleID::gamma ) ++np;
  if ( ng < 2 || ng + np != 3 ) return false;
  return true;
}

ParticleVector Onium3GDecayer::decay(const DecayMode & dm,
				     const Particle & parent) const {
  PVector children = FlatDecayer::decay(dm, parent);
  PVector gluons;
  for ( int i = 0, N = children.size(); i < N; ++i ) {
    children[i]->scale(sqr(parent.mass()));
    if ( children[i]->id() == ParticleID::g ) gluons.push_back(children[i]);
  }
  for ( int i = 0, N = gluons.size(); i < N; ++i )
    gluons[i]->colourNeighbour(gluons[(i + 1)%N]);
  HintPtr h = ptr_new<HintPtr>();
  h->tag(children.begin(), children.end());
  using namespace Group;
  generator()->currentEventHandler()->
    addStep(main, shower()? cascade: hadron, tStepHdlPtr(), h);

  return children;
}

double Onium3GDecayer::reweight(const DecayMode &, const Particle & parent,
				const ParticleVector & ch) const {
  vector<double> x(3);
  Energy2 m2 = parent.momentum().mass2();
  x[0] = 2.0*ch[1]->momentum()*ch[2]->momentum()/m2;
  x[1] = 2.0*ch[2]->momentum()*ch[0]->momentum()/m2;
  x[2] = 2.0*ch[0]->momentum()*ch[1]->momentum()/m2;
  for ( int i = 0; i < 3; ++i )
    if ( ch[i]->id() == ParticleID::gamma &&
	 (1.0 - x[i])*m2 < sqr(minGGMass()) ) return 0.0;
  return 0.5*(sqr((1.0 - x[0])/(x[1]*x[2])) +
	      sqr((1.0 - x[1])/(x[2]*x[0])) +
	      sqr((1.0 - x[2])/(x[0]*x[1])));
}

void Onium3GDecayer::persistentOutput(PersistentOStream & os) const {
  os << doShower << ounit(theMinGGMass,GeV);
}

void Onium3GDecayer::persistentInput(PersistentIStream & is, int) {
  is >> doShower >> iunit(theMinGGMass,GeV);
}

ClassDescription<Onium3GDecayer> Onium3GDecayer::initOnium3GDecayer;
// Definition of the static class description member.

void Onium3GDecayer::Init() {

  static ClassDocumentation<Onium3GDecayer> documentation
    ("This class performs the decay of a spin-1 onium resonance into "
     "three gluons or two gluons and a photon. After the decay the "
     "collision handler is instructed to restart the generation from the "
     "hadronization (or optionally the parton cascade) stage.");

  static Switch<Onium3GDecayer,bool> interfaceShower
    ("Shower",
     "Should the produced gluons be showered or only hadronized?",
     &Onium3GDecayer::doShower, true, true, false);
  static SwitchOption interfaceShowerYes
    (interfaceShower,
     "Yes",
     "The produced gluons should be showered before hadronization.",
     true);
  static SwitchOption interfaceShowerNo
    (interfaceShower,
     "No",
     "The produced gluons should be hadronized whithout preceeding shower.",
     false);

  static Parameter<Onium3GDecayer,Energy> interfaceMinGGMass
    ("MinGGMass",
     "The minimum invariant mass of the two gluons allowed in gamma-g-g "
     "decays.",
     &Onium3GDecayer::theMinGGMass, GeV, 2.0*GeV, ZERO, 10.0*GeV,
     true, false, true);

  interfaceShower.rank(10);

}

