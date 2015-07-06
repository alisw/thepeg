// -*- C++ -*-
//
// Decayer.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Decayer class.
//

#include "Decayer.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace ThePEG;

double Decayer::brat(const DecayMode &,
		     const ParticleData &, double b) const {
  return b;

}

ParticleVector Decayer::getChildren(const DecayMode & dm,
					const Particle &) const {
  return dm.produceProducts();
}

void Decayer::finalBoost(const Particle & parent,
			 const ParticleVector & children) const {
  Utilities::setMomentum(children.begin(), children.end(),
			 parent.momentum().vect(), 1.0e-12);
}  

void Decayer::setScales(const Particle & parent,
			 const ParticleVector & children) const {
  for ( ParticleVector::size_type i = 0; i < children.size(); ++i )
    children[i]->scale(parent.momentum().mass2());
}  

double Decayer::brat(const DecayMode &, const Particle &, double b) const {
  return b;
}

bool Decayer::needsFullStep() const {
  return false;
}

ParticleVector Decayer::
decay(const DecayMode & dm, const Particle & p, Step &) const {
  return decay(dm, p);
}

AbstractNoPIOClassDescription<Decayer> Decayer::initDecayer;

void Decayer::Init() {

  static ClassDocumentation<Decayer> documentation
    ("There is no documentation for the ThePEG::Decayer class");

  static Reference<Decayer,Amplitude> interfaceAmplitude
    ("Amplitude",
     "The eventual amplitude associated to this decay matrix element.",
     &Decayer::theAmplitude, false, false, true, true);

}

void Decayer::persistentOutput(PersistentOStream & os) const {
  os << theAmplitude;
}

void Decayer::persistentInput(PersistentIStream & is, int) {
  is >> theAmplitude;
}

ParticleVector Decayer::DecayParticle(tPPtr parent, Step & s, long maxtry) {
  ParticleVector children;
  if ( !parent ) return children;
  parent = parent->final();
  if ( parent->decayed() ) return children;
  long itry = 0;
  while ( true ) {
    if ( itry++ >=  maxtry ) Throw<DecayFailure>()
      << "Could not decay particle " << parent->data().PDGName() << " after "
      << maxtry << " attempts. Giving up." << Exception::eventerror;
    tDMPtr dm = parent->data().selectMode(*parent);
    if ( !dm ) Throw<DecayFailure>()
      << "Could not decay particle " << parent->data().PDGName() << " since "
      << " no decay mode was found." << Exception::runerror;
    if ( !dm->decayer() ) Throw<DecayFailure>()
      << "Could not perform the decay " << dm->tag()
      << " since the decay mode was not associated with a decayer."
      << Exception::runerror;
    try {
      if ( dm->decayer()->needsFullStep() )
	children = dm->decayer()->decay(*dm, *parent, s);
      else
	children = dm->decayer()->decay(*dm, *parent);
      if ( !children.empty() ) {
	parent->decayMode(dm);
	for ( int i = 0, N = children.size(); i < N; ++i )
	  if ( !s.addDecayProduct(parent, children[i]) )Throw<DecayFailure>()
	    << "An error occurred when tryin to decay an unstable particle "
	    << "of type " << parent->data().PDGName() << ". One of the "
	    << "produced children (of type " << children[i]->data().PDGName()
	    << ") could not be added to the current step."
	    << Exception::abortnow;
	parent->scale(ZERO);
	return children;
      }
    }
    catch (DecayFailure & e) {
      throw e;
    }
    catch (Veto) {}
  }
}
