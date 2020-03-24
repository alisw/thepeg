// -*- C++ -*-
//
// SimpleZGenerator.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleZGenerator class.
//

#include "SimpleZGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

SimpleZGenerator::~SimpleZGenerator() {}

IBPtr SimpleZGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr SimpleZGenerator::fullclone() const {
  return new_ptr(*this);
}

double SimpleZGenerator::generate(cPDPtr q1, cPDPtr q2, Energy2) const {
  if ( BaryonMatcher::Check(*q1) || DiquarkMatcher::Check(*q1) ) {
    if ( BaryonMatcher::Check(*q2) || DiquarkMatcher::Check(*q2) )
      return rnd();
    else
      return sqrt(rnd());
  } else {
    if ( BaryonMatcher::Check(*q2) || DiquarkMatcher::Check(*q2) )
      return 1.0 - sqrt(rnd());
    else
      return rnd();
  }
}

void SimpleZGenerator::persistentOutput(PersistentOStream &) const {}

void SimpleZGenerator::persistentInput(PersistentIStream &, int) {}

ClassDescription<SimpleZGenerator> SimpleZGenerator::initSimpleZGenerator;
// Definition of the static class description member.

void SimpleZGenerator::Init() {

  static ClassDocumentation<SimpleZGenerator> documentation
    ("Implements a naive unphysical model to generate the momentum fraction "
     "\\f$z\\f$ taken by hadrons produced in a hadronization scenario. It "
     "should only be used for testing purposes.");

}

