// -*- C++ -*-
//
// SSSSVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSSSVertex class.
//

#include "SSSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<SSSSVertex> SSSSVertex::initSSSSVertex;
// Definition of the static class description member.

void SSSSVertex::Init() {
  
  static ClassDocumentation<SSSSVertex> documentation
    ("The SSSSVertex class is the implementation"
     "of the helicity amplitude for the four scalar vertex"
     "all vertices of trhis type should inherit from it");
}

// evaluate the vertex
Complex SSSSVertex::evaluate(Energy2 q2, const ScalarWaveFunction & sca1,
			     const ScalarWaveFunction & sca2, 
			     const ScalarWaveFunction & sca3, 
			     const ScalarWaveFunction & sca4) {
  // calculate the coupling
  setCoupling(q2,sca1.particle(),sca2.particle(),
	      sca3.particle(),sca4.particle());
  // return the answer
  return Complex(0.,1.)*norm()*sca1.wave()*sca2.wave()*sca3.wave()*sca4.wave();
}

// off-shell scalar
ScalarWaveFunction SSSSVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
					const ScalarWaveFunction & sca1,
					const ScalarWaveFunction & sca2,
					const ScalarWaveFunction & sca3) {
  // outgoing momentum 
  Lorentz5Momentum pout = sca1.momentum()+sca2.momentum()+sca3.momentum();
  // calculate the coupling
  setCoupling(q2,sca1.particle(),sca2.particle(),sca3.particle(),out);
  // wavefunction
  Energy2 p2   = pout.m2();
  Complex fact = -norm()*sca1.wave()*sca2.wave()*sca3.wave()*
    propagator(iopt,p2,out);
  return ScalarWaveFunction(pout,out,fact);
}
