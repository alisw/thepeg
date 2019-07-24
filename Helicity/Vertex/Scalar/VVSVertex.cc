// -*- C++ -*-
//
// VVSVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVSVertex class.
//

#include "VVSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<VVSVertex,AbstractVVSVertex>
describeThePEGVVSVertex("ThePEG::VVSVertex", "libThePEG.so");
    
void VVSVertex::Init() {
  
  static ClassDocumentation<VVSVertex> documentation
    ("The VVSVertex class is the implementation of the"
     "vector-vector-scalar vertex. All such vertices should inherit"
     "from it.");
  
}

// evaluate the vertex
Complex VVSVertex::evaluate(Energy2 q2,const VectorWaveFunction & vec1,
			    const VectorWaveFunction & vec2, 
			    const ScalarWaveFunction & sca) {
  // calculate the coupling
  setCoupling(q2,vec1.particle(),vec2.particle(),sca.particle());
  // evaluate the vertex
  return Complex(0.,1.)*norm()*sca.wave()*vec1.wave().dot(vec2.wave());
}

// evaluate an off-shell vector
VectorWaveFunction VVSVertex::evaluate(Energy2 q2, int iopt,tcPDPtr out,
				       const VectorWaveFunction & vec,
				       const ScalarWaveFunction & sca,
				       complex<Energy> mass,
				       complex<Energy> width) {
  // outgoing momentum 
  Lorentz5Momentum pout = vec.momentum()+sca.momentum();
  // calculate the coupling
  setCoupling(q2,out,vec.particle(),sca.particle());
  // prefactor
  Energy2 p2    = pout.m2();
  if(mass.real() < ZERO) mass   = out->mass();
  complex<Energy2> mass2 = sqr(mass);
  Complex fact  = norm()*sca.wave()*propagator(iopt,p2,out,mass,width);
  // evaluate the wavefunction
  LorentzPolarizationVector vect;
  // massless case
  if(mass.real()==ZERO) {
    vect = fact*vec.wave();
  }
  // massive case
  else {
    complex<InvEnergy> dot = vec.wave().dot(pout)/mass2;
    vect = fact*(vec.wave()-dot*pout);
  }
  return VectorWaveFunction(pout,out,vect);
}

// off-shell scalar
ScalarWaveFunction VVSVertex::evaluate(Energy2 q2, int iopt,tcPDPtr out, 
				       const VectorWaveFunction & vec1,
				       const VectorWaveFunction & vec2,
				       complex<Energy> mass,
				       complex<Energy> width) {
  // outgoing momentum 
  Lorentz5Momentum pout = vec1.momentum()+vec2.momentum();
  // calculate the coupling
  setCoupling(q2,vec1.particle(),vec2.particle(),out);
  // prefactor
  Energy2 p2   =  pout.m2();
  Complex fact = -norm()*propagator(iopt,p2,out,mass,width);
  // evaluate the wavefunction
  Complex output = fact*vec1.wave().dot(vec2.wave());
  return ScalarWaveFunction(pout,out,output);
}
