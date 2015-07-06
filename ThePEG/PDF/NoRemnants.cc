// -*- C++ -*-
//
// NoRemnants.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NoRemnants class.
//

#include "NoRemnants.h"
#include "ThePEG/Vectors/LorentzVector.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/EventRecord/Particle.h"

using namespace ThePEG;

IBPtr NoRemnants::clone() const {
  return new_ptr(*this);
}

IBPtr NoRemnants::fullclone() const {
  return new_ptr(*this);
}

Lorentz5Momentum NoRemnants::
generate(PartonBinInstance & pb, const double *,
	 Energy2, const LorentzMomentum & p, bool) const {
  if ( pb.particleData() != pb.partonData() || pb.li() > 0.0 )
    throw RemnantHandlerException
      (pb.particleData()->name(), pb.partonData()->name(), name(),
       "This remnant handler cannot handle any remnants!");
  return p;
}

Lorentz5Momentum NoRemnants::
generate(PartonBinInstance & pb, const double *, Energy2,
	 Energy2, const LorentzMomentum & p, bool) const {
  if ( pb.particleData() != pb.partonData() || pb.li() > 0.0 )
    throw RemnantHandlerException
      (pb.particleData()->name(), pb.partonData()->name(), name(),
       "This remnant handler cannot handle any remnants!");
  return p;
}

NoPIOClassDescription<NoRemnants> NoRemnants::initNoRemnants;

void NoRemnants::Init() {}

