// -*- C++ -*-
//
// ColourPairDecayer.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourPairDecayer class.
//

#include "ColourPairDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/HandlerGroup.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/Handlers/EventHandler.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

IBPtr ColourPairDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr ColourPairDecayer::fullclone() const {
  return new_ptr(*this);
}

bool ColourPairDecayer::accept(const DecayMode & dm) const {
  if ( !FlatDecayer::accept(dm) ) return false;
  for ( int i = 0, N = dm.orderedProducts().size(); i < N; ++i ) {
    if ( !dm.orderedProducts()[i]->coloured() ) continue;
    if ( i == N - 1 ) return false;
    if ( dm.orderedProducts()[i]->hasColour() &&
	 !dm.orderedProducts()[i + 1]->hasAntiColour() ) return false;
    if ( dm.orderedProducts()[i]->hasAntiColour() &&
	 !dm.orderedProducts()[i + 1]->hasColour() ) return false;
    ++i;
  }
  return true;
}

ParticleVector ColourPairDecayer::getChildren(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  for ( int i = 0, N = children.size(); i < N; ++i ) {
    children[i]->scale(sqr(parent.mass()));
    if ( !children[i]->coloured() ) continue;
    if ( children[i]->hasColour() )
      children[i]->antiColourNeighbour(children[i + 1]);
    if ( children[i]->hasAntiColour() )
      children[i]->colourNeighbour(children[i + 1]);
    ++i;
  }
  HintPtr h = ptr_new<HintPtr>();
  h->tag(children.begin(), children.end());
  using namespace Group;
  generator()->currentEventHandler()->
    addStep(main, shower()? cascade: hadron, tStepHdlPtr(), h);
  return children;
}


void ColourPairDecayer::persistentOutput(PersistentOStream & os) const {
  os << doShower;
}

void ColourPairDecayer::persistentInput(PersistentIStream & is, int) {
  is >> doShower;
}

ClassDescription<ColourPairDecayer> ColourPairDecayer::initColourPairDecayer;
// Definition of the static class description member.

void ColourPairDecayer::Init() {

  static ClassDocumentation<ColourPairDecayer> documentation
    ("This class performs decays according to phase space into two or "
     "more particles, some of which may be coloured. The coloured "
     "particles must come in pairs and will be colour connected "
     "pair-wise.");

  static Switch<ColourPairDecayer,bool> interfaceShower
    ("Shower",
     "Should the produced partons be showered or only hadronized?",
     &ColourPairDecayer::doShower, true, true, false);
  static SwitchOption interfaceShowerYes
    (interfaceShower,
     "Yes",
     "The produced partons should be showered before hadronization.",
     true);
  static SwitchOption interfaceShowerNo
    (interfaceShower,
     "No",
     "The produced partons should be hadronized whithout preceeding shower.",
     false);

  interfaceShower.rank(10);

}

