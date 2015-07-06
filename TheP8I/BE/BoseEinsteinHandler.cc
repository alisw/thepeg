// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BoseEinsteinHandler class.
//

#include "BoseEinsteinHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/RemnantDecayer.h"   
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/Handlers/EventHandler.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace TheP8I;

BoseEinsteinHandler::BoseEinsteinHandler()
  :  pythia(), doContinueDecays(true),
#include "BoseEinsteinHandler-init.h"
      {}

BoseEinsteinHandler::~BoseEinsteinHandler() {}

void BoseEinsteinHandler::
performDecay(tPPtr parent, Step & s) const {
  if ( parent->data().width() < theBoseEinstein_widthSep*GeV &&
       parent->id() != ParticleID::K0 &&
       parent->id() != ParticleID::Kbar0 ) return;
  ParticleVector children = Decayer::DecayParticle(parent, s, maxLoop());
  for ( int i = 0, N = children.size(); i < N; ++i )
    if ( !children[i]->data().stable() ) performDecay(children[i], s);
}

void BoseEinsteinHandler::handle(EventHandler & eh, const tPVector & tagged,
	    const Hint & h) {
  if ( !pythia.created() ) {
    pythia.init(*this, theAdditionalP8Settings);
    theBE.init(&(pythia().info), pythia().settings, pythia().particleData);
  }

  // First go through to see which of the tagged particles should be
  // decayed and/or shifted: Exit if none are found.
  tPVector parents;
  for ( int i = 0, N = tagged.size(); i < N; ++i )
    if ( tagged[i] ) parents.push_back(tagged[i]);

  if ( parents.empty() ) return;

  // Create a new step, decay all particles which have a width larger
  // than the given limit
  // (<interface>BoseEinstein_widthSep</interface>) and add their
  // children in the new step.
  for ( int i = 0, N = parents.size(); i < N; ++i )
    if ( !parents[i]->data().stable() )
      performDecay(newStep()->find(parents[i]->final()), *newStep());

  // Put all particles (or their children fi they have decayed) in the
  // Pythia8::Event.
  pythia.clearEvent();
  particles.clear();
  for ( int i = 0, N = parents.size(); i < N; ++i ) addParticle(parents[i]);

  // Perform the BE-shift
  theBE.shiftEvent(pythia.event());

  // Insert all shifted particles in the Step.
  for ( int i = 0, N = pythia.event().size(); i < N; ++i ) {
    if ( !pythia.event()[i].isFinal() ) continue;
    int im = pythia.event()[i].mother1();
    if ( im <= 0 ) continue;
    Lorentz5Momentum p(pythia.event()[i].px()*GeV,
		       pythia.event()[i].py()*GeV,
		       pythia.event()[i].pz()*GeV,
		       pythia.event()[i].e()*GeV,
		       pythia.event()[i].m()*GeV);
    if ( particles[im]->momentum().x() != p.x() ||
	 particles[im]->momentum().y() != p.y() ||
	 particles[im]->momentum().z() != p.z() ||
	 particles[im]->momentum().e() != p.e() ) {
      tPPtr newp = newStep()->copyParticle(particles[im]);
      newp->set5Momentum(p);
      particles[im] = newp;
    }
  }

  // Now go through all particles and decay them further if requested.
  if ( !doContinueDecays ) return;
  for ( int i = 0, N = particles.size(); i < N; ++i )
    if ( particles[i] && !particles[i]->data().stable() )
      DecayHandler::performDecay(particles[i], *newStep());

}

void BoseEinsteinHandler::addParticle(tPPtr p) {
  if ( p->next() ) addParticle(p->next());
  else if ( p->decayed() )
    for ( int i = 0, N = p->children().size(); i < N; ++i )
      addParticle(p->children()[i]);
  else {
    int idx = pythia.addParticle(p, 1, 0, 0);
    particles.resize(idx + 1);
    particles[idx] = p;
  }
}

IBPtr BoseEinsteinHandler::clone() const {
  return new_ptr(*this);
}

IBPtr BoseEinsteinHandler::fullclone() const {
  return new_ptr(*this);
}


void BoseEinsteinHandler::doinitrun() {
  DecayHandler::doinitrun();
  theAdditionalP8Settings.push_back("ProcessLevel:all = off");
  theAdditionalP8Settings.push_back("HadronLevel:Decay = off");
  theAdditionalP8Settings.push_back("Check:event = off");
  pythia.init(*this, theAdditionalP8Settings);
  theBE.init(&(pythia().info), pythia().settings, pythia().particleData);
}

void BoseEinsteinHandler::persistentOutput(PersistentOStream & os) const {
  os
#include "BoseEinsteinHandler-output.h"
    << doContinueDecays;
}

void BoseEinsteinHandler::persistentInput(PersistentIStream & is, int) {
  is
#include "BoseEinsteinHandler-input.h"
    >> doContinueDecays;
}

ClassDescription<BoseEinsteinHandler> BoseEinsteinHandler::initBoseEinsteinHandler;
// Definition of the static class description member.

void BoseEinsteinHandler::Init() {

#include "BoseEinsteinHandler-interfaces.h"



  static Switch<BoseEinsteinHandler,bool> interfaceContinueDecays
    ("ContinueDecays",
     "Continue the decay chains after the Bose-Einstein shifting.",
     &BoseEinsteinHandler::doContinueDecays, true, true, false);
  static SwitchOption interfaceContinueDecaysYes
    (interfaceContinueDecays,
     "Yes",
     "Continue decays.",
     true);
  static SwitchOption interfaceContinueDecaysNo
    (interfaceContinueDecays,
     "No",
     "Do not continue decays.",
     false);


}

