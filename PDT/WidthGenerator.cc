// -*- C++ -*-
//
// WidthGenerator.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WidthGenerator class.
//

#include "WidthGenerator.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/UseRandom.h"

using namespace ThePEG;

Length WidthGenerator::lifeTime(const ParticleData &, Energy, Energy w) const {
  return UseRandom::rndExp(hbarc/w);
}

WidthGenerator::DecayMap WidthGenerator::rate(const Particle & p) {
  DecayMap dm(p.data().decaySelector());
  return dm;
}

Energy WidthGenerator::width(const Particle & p) const {
  return width(p.data(), p.mass());
}

AbstractNoPIOClassDescription<WidthGenerator>
WidthGenerator::initWidthGenerator;

void WidthGenerator::Init() {

  static ClassDocumentation<WidthGenerator> documentation
    ("Abstract base class to be used to encapsulate models for giving the "
     "partial decay width of a ParticleData given the specified DecayModes.");

}

