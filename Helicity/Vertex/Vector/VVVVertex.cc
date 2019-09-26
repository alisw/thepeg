// -*- C++ -*-
//
// VVVVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVVVertex class.
//

#include "VVVVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity; 

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<VVVVertex,AbstractVVVVertex>
describeThePEGVVVVertex("ThePEG::VVVVertex", "libThePEG.so");
  
void VVVVertex::Init() {
  
  static ClassDocumentation<VVVVertex> documentation
    ("The VVVVertex class implements the helicity amplitude"
     "calculations for the triple gauge boson vertex. Any   "
     "implementation of such a vertex should inherit from in and implement"
     " the virtual setCoupling member to calculate the coupling");
  
}

// evaluate the vertex
Complex VVVVertex::evaluate(Energy2 q2, const VectorWaveFunction & vec1,
			    const VectorWaveFunction & vec2,
			    const VectorWaveFunction & vec3) {
  // calculate the coupling
  setCoupling(q2,vec1.particle(),vec2.particle(),vec3.particle());
  complex<Energy> alpha1(ZERO);
  // decide if we need to use special treatment to avoid gauge cancelations
  // first vector
  if(abs(vec1.t())!=0.) {
    if(abs(vec1.t())>0.1*max( max(abs(vec1.x()),abs(vec1.y())),abs(vec1.z())))
      alpha1=vec1.e()/vec1.t();
  }
  // second vector
  if(abs(vec2.t())!=0.) {
    if(abs(vec2.t())>0.1*max( max(abs(vec2.x()),abs(vec2.y())),abs(vec2.z())))
      alpha1=vec2.e()/vec2.t();
  }
  // third vector
  if(abs(vec3.t())!=0.) {
    if(abs(vec3.t())>0.1*max( max(abs(vec3.x()),abs(vec3.y())),abs(vec3.z())))
      alpha1=vec3.e()/vec3.t();
  }
  // dot products of the polarization vectors
  Complex dot12 = vec1.wave().dot(vec2.wave());
  Complex dot13 = vec1.wave().dot(vec3.wave());
  Complex dot23 = vec3.wave().dot(vec2.wave());
  // dot products of polarization vectors and momentum
  complex<Energy> dotp13 = vec3.wave().dot(LorentzPolarizationVectorE(vec1.momentum())
					   - alpha1 * vec1.wave());
  complex<Energy> dotp23 = vec3.wave().dot(LorentzPolarizationVectorE(vec2.momentum())
					   - alpha1 * vec2.wave());
  complex<Energy> dotp21 = vec1.wave().dot(LorentzPolarizationVectorE(vec2.momentum())
					   - alpha1 * vec2.wave());
  complex<Energy> dotp31 = vec1.wave().dot(LorentzPolarizationVectorE(vec3.momentum())
					   - alpha1 * vec3.wave());
  complex<Energy> dotp32 = vec2.wave().dot(LorentzPolarizationVectorE(vec3.momentum())
					   - alpha1 * vec3.wave());
  complex<Energy> dotp12 = vec2.wave().dot(LorentzPolarizationVectorE(vec1.momentum())
					   - alpha1 * vec1.wave());
  // finally calculate the vertex
  return Complex(0.,1.)*norm()*UnitRemoval::InvE*
    (dot12*(dotp13-dotp23)+dot23*(dotp21-dotp31)+dot13*(dotp32-dotp12));
}
  
// off-shell vector
VectorWaveFunction VVVVertex::evaluate(Energy2 q2,int iopt, tcPDPtr out,
				       const VectorWaveFunction & vec1,
				       const VectorWaveFunction & vec2,
				       complex<Energy> mass, complex<Energy> width) {
  // output momenta
  Lorentz5Momentum pout =vec1.momentum()+vec2.momentum();
  // calculate the coupling
  setCoupling(q2,out,vec1.particle(),vec2.particle());
  // prefactor
  Energy2 p2    = pout.m2();
  Complex fact  = norm()*propagator(iopt,p2,out,mass,width);
  if(mass.real() < ZERO) mass   = out->mass();
  complex<Energy2> mass2 = sqr(mass);
  // dot products we need
  Complex dot12 = vec1.wave().dot(vec2.wave());
  complex<Energy> dota = vec1.wave().dot(pout+vec2.momentum());
  complex<Energy> dotb = vec2.wave().dot(pout+vec1.momentum());
  // compute the polarization vector
  LorentzPolarizationVector vect = UnitRemoval::InvE*fact*
    (dot12*(vec1.momentum()-vec2.momentum())-dotb*vec1.wave()+dota*vec2.wave());
  // scalar piece for massive case
  if(mass.real()!=ZERO) {
    complex<InvEnergy> dot = vect.dot(pout)/mass2;
    vect -= dot*pout;       
  }
  return VectorWaveFunction(pout,out,vect);
}
