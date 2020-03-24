// -*- C++ -*-
//
// RemnantParticle.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RemnantParticle class.
//

#include "RemnantParticle.h"
#include "ThePEG/PDT/RemnantData.h"
#include "ThePEG/PDT/RemnantDecayer.h"
#include "ThePEG/EventRecord/MultiColour.h"
#include "ThePEG/Config/algorithm.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

RemnantParticle::
RemnantParticle(const Particle & particle, RemDecPtr decayer, tPPtr parton)
  : Particle(new_ptr(RemnantData(particle.dataPtr(), decayer))) {
  remData = const_ptr_cast<tRemPDPtr>(dynamic_ptr_cast<tcRemPDPtr>(dataPtr()));
  set5Momentum(particle.momentum());
  colourInfo(new_ptr(MultiColour()));
  parent = &particle;
  if ( parton ) extract(parton);
}

bool RemnantParticle::extract(tPPtr parton, bool fixcolour) {
  LorentzMomentum pnew = momentum() - parton->momentum();
  if ( !remData->decayer().checkExtract(parent, parton, pnew) ) return false;
  if ( !remData->extract(parton->dataPtr()) ) return false;
  theExtracted.push_back(parton);
  setMomentum(pnew);
  rescaleMass();
  if ( fixcolour ) fixColourLines(parton);
  return true;
}

bool RemnantParticle::reextract(tPPtr oldp, tPPtr newp, bool fixcolour) {
  LorentzMomentum pnew = momentum() + oldp->momentum() - newp->momentum();
  if ( !remData->decayer().checkExtract(parent, newp, pnew) ) return false;
  PVector::iterator it = find(theExtracted, oldp);
  if ( it == theExtracted.end() ) return false;
  if ( !remData->reextract(oldp->dataPtr(), newp->dataPtr()) ) return false;
  theExtracted[it - theExtracted.begin()] = newp;
  setMomentum(pnew);
  rescaleMass();
  if ( oldp->colourLine() ) oldp->colourLine()->removeAntiColoured(this);
  if ( oldp->antiColourLine() ) oldp->antiColourLine()->removeColoured(this);
  if ( fixcolour )  fixColourLines(newp);
  return true;
}

bool RemnantParticle::remove(tPPtr oldp) {
  LorentzMomentum pnew = momentum() + oldp->momentum();
  PVector::iterator it = find(theExtracted, oldp);
  if ( it == theExtracted.end() ) return false;
  if ( !remData->remove(oldp->dataPtr()) ) return false;
  theExtracted.erase(it);
  setMomentum(pnew);
  rescaleMass();
  if ( oldp->colourLine() ) oldp->colourLine()->removeAntiColoured(this);
  if ( oldp->antiColourLine() ) oldp->antiColourLine()->removeColoured(this);
  return true;
}

void RemnantParticle::fixColourLines(tPPtr parton) {
  if ( parton->hasColour() ) {
    if ( parton->colourLine() )
      parton->colourLine()->addAntiColoured(this);
    else
      ColourLine::create(parton, tPPtr(this));
  } 
  if ( parton->hasAntiColour() ) {
    if ( parton->antiColourLine() )
      parton->antiColourLine()->addColoured(this);
    else ColourLine::create(this, parton);
  }
}

void RemnantParticle::persistentOutput(PersistentOStream & os) const {
  os << remData << parent << theExtracted;
}

void RemnantParticle::persistentInput(PersistentIStream & is, int) {
  is >> remData >> parent >> theExtracted;
}

ClassDescription<RemnantParticle> RemnantParticle::initRemnantParticle;
// Definition of the static class description member.

void RemnantParticle::Init() {}

