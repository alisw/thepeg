// -*- C++ -*-
//
// GaussianPtGenerator.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GaussianPtGenerator class.
//

#include "GaussianPtGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Vectors/Transverse.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

GaussianPtGenerator::~GaussianPtGenerator() {}

IBPtr GaussianPtGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr GaussianPtGenerator::fullclone() const {
  return new_ptr(*this);
}

TransverseMomentum GaussianPtGenerator::generate() const {
  pair<Energy,Energy> ret;
  Energy pt = ZERO;
  while ( ( pt = theSigma*sqrt(-log(rnd())) ) > theUpperCut ) {}
  double phi = rnd(2.0*Constants::pi);
  ret.first = pt*cos(phi);
  ret.second = pt*sin(phi);

  return ret;
}


void GaussianPtGenerator::persistentOutput(PersistentOStream & os) const {
  os << ounit(theSigma, GeV) << ounit(theUpperCut, GeV);
}

void GaussianPtGenerator::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theSigma, GeV) >> iunit(theUpperCut, GeV);
}

ClassDescription<GaussianPtGenerator> GaussianPtGenerator::initGaussianPtGenerator;
// Definition of the static class description member.

void GaussianPtGenerator::Init() {

  static ClassDocumentation<GaussianPtGenerator> documentation
    ("The ThePEG::GaussianPtGenerator class generates a gaussian "
     "transverse momentum.");

  static Parameter<GaussianPtGenerator,Energy> interfaceSigma
    ("Sigma",
     "The width of the Gaussian distribution. The average squared transverse "
     "momentum is Sigma squared.",
     &GaussianPtGenerator::theSigma, GeV, 1.0*GeV, ZERO, 10.0*GeV,
     true, false, true);

  static Parameter<GaussianPtGenerator,Energy> interfaceUpperCut
    ("UpperCut",
     "Upper cutoff for the transverse momentum distribution.",
     &GaussianPtGenerator::theUpperCut, GeV, 5.0*GeV, ZERO, 50.0*GeV,
     true, false, true);

  interfaceSigma.rank(10);

}

