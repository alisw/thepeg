// -*- C++ -*-
//
// DummyDecayer.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DummyDecayer class.
//

#include "DummyDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

IBPtr DummyDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr DummyDecayer::fullclone() const {
  return new_ptr(*this);
}

bool DummyDecayer::accept(const DecayMode &) const {
  return true;
}

ParticleVector DummyDecayer::decay(const DecayMode &,
				  const Particle &) const {
  throw std::logic_error("Tried to decay with the DummyDecayer class.");
}

double DummyDecayer::
brat(const DecayMode &, const ParticleData &, double) const {
  return 0.0;
}

double DummyDecayer::brat(const DecayMode &, const Particle &, double) const {
  return 0.0;
}


NoPIOClassDescription<DummyDecayer> DummyDecayer::initDummyDecayer;
// Definition of the static class description member.

void DummyDecayer::Init() {

  static ClassDocumentation<DummyDecayer> documentation
    ("This is a dummy decayer class to be used for symbolic decay "
     "channels. If it for some reason is called to perform a decay, it "
     "will throw a std::logic_error.");

}

