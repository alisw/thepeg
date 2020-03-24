// -*- C++ -*-
//
// Hint.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Hint class.
//

#include "Hint.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/Particle.h"

using namespace ThePEG;

tPVector Hint::tagged(const Step & s) const {
  tPVector tags;
  if ( stop() ) return tags;
  if ( !tagged() ) {
    tags.insert(tags.end(), s.particles().begin(), s.particles().end());
    return tags;
  }
  for ( int i = 0, N = theTagged.size(); i < N; ++i ) {
    tPPtr p = theTagged[i];
    while ( !s.find(p) && p->next() ) p = p->next();
    if ( s.find(p) ) tags.push_back(p);
  }
  return tags;
}

void Hint::persistentOutput(PersistentOStream & os) const {
  os << ounit(theScale, GeV2) << theStopFlag << theTagged;
}

void Hint::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theScale, GeV2) >> theStopFlag >> theTagged;
}

ClassDescription<Hint> Hint::initHint;

void Hint::Init() {}

Hint Hint::theDefaultHint;
