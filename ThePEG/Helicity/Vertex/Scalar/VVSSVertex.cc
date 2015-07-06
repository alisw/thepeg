// -*- C++ -*-
//
// VVSSVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVSSVertex class.
//

#include "VVSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<VVSSVertex> VVSSVertex::initVVSSVertex;
// Definition of the static class description member.
    
void VVSSVertex::Init() {
      
static ClassDocumentation<VVSSVertex> documentation
  ("The VVSSVertex class is the implementation of helicity"
   "amplitude calculation of the vector-vector-scalar-scalar vertex."
   "All classes for this type of vertex should inherit from it.");
}
 
// evaluate the vertex
Complex VVSSVertex::evaluate(Energy2 q2,const VectorWaveFunction & vec1,
			     const VectorWaveFunction & vec2, 
			     const ScalarWaveFunction & sca1, 
			     const ScalarWaveFunction & sca2) {
  // calculate the coupling
  setCoupling(q2,vec1.particle(),vec2.particle(),
	      sca1.particle(),sca2.particle());
  // evaluate the vertex
  return Complex(0.,1.)*norm()*sca1.wave()*sca2.wave()*
    vec1.wave().dot(vec2.wave());
}

// evaluate an off-shell vector
VectorWaveFunction VVSSVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
					const VectorWaveFunction & vec,
					const ScalarWaveFunction & sca1,
					const ScalarWaveFunction & sca2,
					complex<Energy> mass,
					complex<Energy> width) {
  // outgoing momentum 
  Lorentz5Momentum pout = vec.momentum()+sca1.momentum()+sca2.momentum();
  // calculate the coupling
  setCoupling(q2,out,vec.particle(),sca1.particle(),sca2.particle());
  // prefactor
  Energy2 p2    = pout.m2();
  if(mass.real() < ZERO) mass   = out->mass();
  complex<Energy2> mass2 = sqr(mass);
  Complex fact  = norm()*sca1.wave()*sca2.wave()*propagator(iopt,p2,out,mass,width);
  // evaluate the wavefunction
  LorentzPolarizationVector vect;
  // massless case
  if(mass.real()==ZERO) {
    vect = fact*vec.wave();
  }
  // massive case
  else {
    complex<InvEnergy> dot = vec.wave().dot(pout)/mass2;
    vect = fact*(vec.wave()-(dot*pout));
  }
  return VectorWaveFunction(pout,out,vect);
}

// off-shell scalar
ScalarWaveFunction VVSSVertex::evaluate(Energy2 q2, int iopt,tcPDPtr out, 
					const VectorWaveFunction & vec1,
					const VectorWaveFunction & vec2,
					const ScalarWaveFunction & sca,
					complex<Energy> mass,
					complex<Energy> width) {
  // outgoing momentum 
  Lorentz5Momentum pout = vec1.momentum()+vec2.momentum()+sca.momentum(); 
  // calculate the coupling
  setCoupling(q2,vec1.particle(),vec2.particle(),out,sca.particle());
  // prefactor
  Energy2 p2   =  pout.m2();
  Complex fact = -norm()*sca.wave()*propagator(iopt,p2,out,mass,width);
  // evaluate the wavefunction
  Complex output = fact*vec1.wave().dot(vec2.wave());
  return ScalarWaveFunction(pout,out,output);
}
