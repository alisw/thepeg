// -*- C++ -*-
//
// PDT.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PDT class.
//

#include "PDT.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/EventRecord/Particle.h"

using namespace ThePEG;

vector<long> PDT::flavourContent(long id) {
  vector<long> ret;
  if ( id == ParticleID::K_L0 || id == ParticleID::K_S0 ) {
    ret.push_back(ParticleID::s);
    ret.push_back(ParticleID::d);
  }
  else if ( MesonMatcher::Check(id) ) {
    ret.push_back((id/100)%10);
    ret.push_back(-(id/10)%10);
  }
  else if ( BaryonMatcher::Check(id) ) {
    ret.push_back((id/1000)%10);
    long iqb = (id/100)%10;
    long iqc = (id/10)%10;
    if ( abs(iqb) < abs(iqc) ) {
      ret.push_back(iqc);
      ret.push_back(iqb);
    } else {
      ret.push_back(iqb);
      ret.push_back(iqc);
    }
  }
  else if ( DiquarkMatcher::Check(id) ) {
    ret.push_back((id/1000)%10);
    ret.push_back((id/100)%10);
  }
  else if ( abs(id) < 10 ) {
    ret.push_back(id);
  }
  return ret;
}

vector<long> PDT::flavourContent(tcPDPtr p) {
  return flavourContent(p->id());
}

vector<long> PDT::flavourContent(tcPPtr p) {
  return flavourContent(p->id());
}

vector<long> PDT::flavourContent(const ParticleData & p) {
  return flavourContent(p.id());
}

vector<long> PDT::flavourContent(const Particle & p) {
  return flavourContent(p.id());
}


