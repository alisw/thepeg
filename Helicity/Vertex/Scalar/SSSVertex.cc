// -*- C++ -*-
//
// SSSVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSSVertex class.
//

#include "SSSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<SSSVertex,AbstractSSSVertex>
describeThePEGSSSVertex("ThePEG::SSSVertex", "libThePEG.so");

void SSSVertex::Init() {
  
  static ClassDocumentation<SSSVertex> documentation
    ("The SSSVertex class is the implementation of the SSS"
     "vertex. All such vertices shoud inherit from it");
}

// evaluate the vertex
Complex SSSVertex::evaluate(Energy2 q2,
			    const ScalarWaveFunction & sca1,
			    const ScalarWaveFunction & sca2,
 			    const ScalarWaveFunction & sca3) {
  if(kinematics()) calculateKinematics(sca1.momentum(),sca2.momentum(),sca3.momentum());
  // calculate the coupling
  setCoupling(q2,sca1.particle(),sca2.particle(),sca3.particle());
  // return the answer
  return Complex(0.,1.)*norm()*sca1.wave()*sca2.wave()*sca3.wave();
}

// off-shell scalar
ScalarWaveFunction SSSVertex::evaluate(Energy2 q2,int iopt, tcPDPtr out, 
				       const ScalarWaveFunction & sca1,
				       const ScalarWaveFunction & sca2,
				       complex<Energy> mass,
				       complex<Energy> width) {
  // outgoing momentum 
  Lorentz5Momentum pout = sca1.momentum()+sca2.momentum();
  if(kinematics()) calculateKinematics(sca1.momentum(),sca2.momentum(),-pout);
  // calculate the coupling
  setCoupling(q2,sca1.particle(),sca2.particle(),out);
  // wavefunction
  Energy2 p2=pout.m2();
  Complex fact=-norm()*sca1.wave()*sca2.wave()*propagator(iopt,p2,out,mass,width);
  return ScalarWaveFunction(pout,out,fact);
}
