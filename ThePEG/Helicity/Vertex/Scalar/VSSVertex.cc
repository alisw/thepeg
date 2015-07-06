// -*- C++ -*-
//
// VSSVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VSSVertex class.
//

#include "VSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;
using namespace Helicity;
    
AbstractNoPIOClassDescription<VSSVertex> VSSVertex::initVSSVertex;
// Definition of the static class description member.

void VSSVertex::Init() {
      
static ClassDocumentation<VSSVertex> documentation
  ("The VSSVertex class is hte implementation of the"
   "vector-scalar-scalar vertex for helicity amplitude calculations."
   " all such vertices should inherit from it");
 
}

// evaluate the vertex
Complex VSSVertex::evaluate(Energy2 q2, const VectorWaveFunction & vec,
			    const ScalarWaveFunction & sca1,
			    const ScalarWaveFunction & sca2) {
  // calculate the coupling
  setCoupling(q2,vec.particle(),sca1.particle(),sca2.particle());
  // calculate the vertex
  return UnitRemoval::InvE * -Complex(0.,1.) * norm() * sca1.wave()*sca2.wave()*
    vec.wave().dot(sca1.momentum()-sca2.momentum());
}

// off-shell vector
VectorWaveFunction VSSVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
				       const ScalarWaveFunction & sca1,
				       const ScalarWaveFunction & sca2,
				       complex<Energy> mass,
				       complex<Energy> width) {
  // outgoing momentum 
  Lorentz5Momentum pout(sca1.momentum()+sca2.momentum());
  // calculate the coupling
  setCoupling(q2,out,sca1.particle(),sca2.particle());
  // mass and width
  if(mass.real() < ZERO)  mass   = out->mass();
  complex<Energy2> mass2 = sqr(mass);
  // calculate the prefactor
  Energy2 p2    = pout.m2();
  Complex fact = norm()*sca1.wave()*sca2.wave()*propagator(iopt,p2,out,mass,width);
  // compute the vector
  LorentzPolarizationVector vec;
  // massless outgoing vector
  if(mass.real()==ZERO) {
    vec = UnitRemoval::InvE * fact * (sca2.momentum()-sca1.momentum());
  }
  // massive outgoing vector
  else {
    // first the dot product for the second term
    Complex dot = (sca1.m2()-sca2.m2())/mass2;
    // compute the vector
    vec = fact * 
      (LorentzPolarizationVector(UnitRemoval::InvE * (sca2.momentum()-sca1.momentum()))
       +dot*UnitRemoval::InvE * pout);
  }
  return VectorWaveFunction(pout,out,vec);
}

// return an off-shell scalar
ScalarWaveFunction VSSVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
				       const VectorWaveFunction & vec,
				       const ScalarWaveFunction & sca,
				       complex<Energy> mass,
				       complex<Energy> width ) {
  // momentum of the particle
  Lorentz5Momentum pout = sca.momentum()+vec.momentum(); 
  // calculate the coupling
  setCoupling(q2,vec.particle(),sca.particle(),out);
  // calculate the prefactor
  Energy2 p2   = pout.m2();
  Complex fact = norm()*sca.wave()*propagator(iopt,p2,out,mass,width);
  // compute the wavefunction
  fact = UnitRemoval::InvE * fact*vec.wave().dot(sca.momentum()+pout);
  return ScalarWaveFunction(pout,out,fact);
}
