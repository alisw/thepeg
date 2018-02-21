// -*- C++ -*-
//
// BreitWignerMass.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BreitWignerMass class.
//

#include "BreitWignerMass.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/UseRandom.h"

using namespace ThePEG;

IBPtr BreitWignerMass::clone() const {
  return new_ptr(*this);
}

IBPtr BreitWignerMass::fullclone() const {
  return new_ptr(*this);
}

Energy BreitWignerMass::mass(const ParticleData & pd) const {
  Energy ret = ZERO;
  do {
    ret = UseRandom::rndRelBW(pd.mass(), pd.width(), pd.widthCut());
  } while ( ret > pd.massMax() || ret < pd.massMin() );
  return ret;
}

NoPIOClassDescription<BreitWignerMass> BreitWignerMass::initBreitWignerMass;

void BreitWignerMass::Init() {

  static ClassDocumentation<BreitWignerMass> documentation
    ("Generates masses of particle instances according to a Breit-Wigner "
     "distribution.");

}

