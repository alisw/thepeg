// -*- C++ -*-
//
// FlatDecayer.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FlatDecayer class.
//

#include "FlatDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Vectors/LorentzRotation.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

IBPtr FlatDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FlatDecayer::fullclone() const {
  return new_ptr(*this);
}

bool FlatDecayer::accept(const DecayMode & dm) const {
  if ( dm.products().size() == 1 && 
       ( dm.parent()->massMax() > (**(dm.products().begin())).massMax() ||
	 dm.parent()->massMin() < (**(dm.products().begin())).massMin() ) )
    return false;
  return dm.products().size() > 0 && dm.cascadeProducts().empty() &&
    dm.productMatchers().empty() && !dm.wildProductMatcher();
}

ParticleVector FlatDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector children = getChildren(dm, parent);
  try {
    do {
      if ( children.size() == 1 ) {
	children[0]->setMomentum(parent.momentum());
	children[0]->scale(parent.momentum().mass2());
	return children;
      }
      else {
	SimplePhaseSpace::CMSn(children, parent.mass());
      }
    } while ( reweight(dm, parent, children) < UseRandom::rnd() );
  }
  catch ( ImpossibleKinematics & ) {
    children.clear();
    return children;
  }

  finalBoost(parent, children);
  setScales(parent, children);

  return children;
}

NoPIOClassDescription<FlatDecayer> FlatDecayer::initFlatDecayer;
// Definition of the static class description member.

void FlatDecayer::Init() {

  static ClassDocumentation<FlatDecayer> documentation
    ("The ThePEG::FlatDecayer class describes the decay of a "
     "ThePEG::Particle into a set of specified children according "
     "to a flat distribution in phase space.");

}

