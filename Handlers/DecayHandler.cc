// -*- C++ -*-
//
// DecayHandler.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayHandler class.
//

#include "DecayHandler.h"
#include "Hint.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/Decayer.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

DecayHandler::~DecayHandler() {}

IBPtr DecayHandler::clone() const {
  return new_ptr(*this);
}

IBPtr DecayHandler::fullclone() const {
  return new_ptr(*this);
}

void DecayHandler::
handle(EventHandler &, const tPVector & tagged,
       const Hint &) {
  // First go through to see which of the tagged particles should be
  // decayed: Exit if none are found.
  tPVector parents;
  for ( int i = 0, N = tagged.size(); i < N; ++i )
    if ( tagged[i] && !tagged[i]->data().stable() )
      parents.push_back(tagged[i]);

  if ( parents.empty() ) return;

  // Create a new step, decay all particles and add their children in
  // the new step.
  for ( int i = 0, N = parents.size(); i < N; ++i )
    performDecay(newStep()->find(parents[i]->final()), *newStep());
}

void DecayHandler::
performDecay(tPPtr parent, Step & s) const {
  if ( maxLifeTime() >= ZERO ) {
    if( ( lifeTimeOption() && parent->lifeLength().tau() > maxLifeTime())||
	(!lifeTimeOption() && parent->data().cTau()      > maxLifeTime()) ) {
      parent->setLifeLength(Distance());
      return;
    }
  }
  ParticleVector children = Decayer::DecayParticle(parent, s, maxLoop());
  for ( int i = 0, N = children.size(); i < N; ++i )
    if ( !children[i]->data().stable() ) performDecay(children[i], s);
}

void DecayHandler::persistentOutput(PersistentOStream & os) const {
  os << theMaxLoop << ounit(theMaxLifeTime, mm) << theLifeTimeOption;
}

void DecayHandler::persistentInput(PersistentIStream & is, int) {
  is >> theMaxLoop >> iunit(theMaxLifeTime, mm) >> theLifeTimeOption;
}

ClassDescription<DecayHandler> DecayHandler::initDecayHandler;

void DecayHandler::Init() {

  static ClassDocumentation<DecayHandler> documentation
    ("This is the main class handling the decay of unstable particles. Note "
     "that the actual decays will be performed by objects of the "
     "ThePEG::Decayer class.");

  static Parameter<DecayHandler,long> interfaceMaxLoop
    ("MaxLoop",
     "The maximum number of attempts per event when selecting a decay channel.",
     &DecayHandler::theMaxLoop, 100000, 100, 100000000, true, false, true);

  static Parameter<DecayHandler,Length> interfaceMaxLifeTime
    ("MaxLifeTime",
     "The maximum lifetime (c*tau) in units of mm. Particles with longer "
     "lifetimes than this will not be decayed.",
     &DecayHandler::theMaxLifeTime, mm, -1.0*mm, -1.0*mm, ZERO,
     true, false, Interface::lowerlim);


  static Switch<DecayHandler,bool> interfaceLifeTimeOption
    ("LifeTimeOption",
     "Option for how the maximum life time is interpreted",
     &DecayHandler::theLifeTimeOption, false, false, false);
  static SwitchOption interfaceLifeTimeOptionAverage
    (interfaceLifeTimeOption,
     "Average",
     "Cut on the average lifetime of the particle type",
     false);
  static SwitchOption interfaceLifeTimeOptionReal
    (interfaceLifeTimeOption,
     "Real",
     "Cut on the lifetime generated for the given instance",
     true);

}

