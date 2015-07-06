// -*- C++ -*-
//
// FFSVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFSVertex class.
//

#include "FFSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

AbstractNoPIOClassDescription<FFSVertex> FFSVertex::initFFSVertex;
// Definition of the static class description member.
    
void FFSVertex::Init() {

  static ClassDocumentation<FFSVertex> documentation
    ("The FFSVertex class is the implementation of the FFS"
     "vertex. All such vertices shoud inherit from it");
  
}
  
// evaluate the full vertex
Complex FFSVertex::evaluate(Energy2 q2, const SpinorWaveFunction & sp,
			    const SpinorBarWaveFunction & sbar,
			    const ScalarWaveFunction & sca) {
  // calculate the couplings
  setCoupling(q2,sp.particle(),sbar.particle(),sca.particle());
  Complex vertex(  _left*(sbar.s1()*sp.s1()+sbar.s2()*sp.s2())
		   +_right*(sbar.s3()*sp.s3()+sbar.s4()*sp.s4())
		   );
  // final factors
  return Complex(0.,1.)*norm()*sca.wave()*vertex;
}

// off-shell scalar
ScalarWaveFunction FFSVertex::evaluate(Energy2 q2,int iopt, tcPDPtr out, 
				       const SpinorWaveFunction & sp,
				       const SpinorBarWaveFunction & sbar,
				       complex<Energy> mass,
				       complex<Energy> width) {
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = sbar.momentum()+sp.momentum();
  // first calculate the couplings
  setCoupling(q2,sp.particle(),sbar.particle(),out);
  Energy2 p2   = pout.m2();
  Complex fact = -norm()*propagator(iopt,p2,out,mass,width);
  Complex output =  _left*(sbar.s1()*sp.s1()+sbar.s2()*sp.s2())
    +_right*(sbar.s3()*sp.s3()+sbar.s4()*sp.s4());
  // final factors and output
  output*=fact;
  return ScalarWaveFunction(pout,out,output);
}
    
// off-shell spinor
SpinorWaveFunction FFSVertex::evaluate(Energy2 q2, int iopt,tcPDPtr out,
				       const SpinorWaveFunction & sp,
				       const ScalarWaveFunction & sca,
				       complex<Energy> mass,
				       complex<Energy> width) {
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = sp.momentum()+sca.momentum();
  // first calculate the couplings
  setCoupling(q2,sp.particle(),out,sca.particle());
  Energy2 p2   = pout.m2();
  Complex fact = -norm()*sca.wave()*propagator(iopt,p2,out,mass,width);
  Complex ii(0.,1.);
  // useful combinations of the momenta
  if(mass.real() < ZERO) mass  = out->mass();
  complex<Energy> p1p2 = pout.x()+ii*pout.y();
  complex<Energy> p1m2 = pout.x()-ii*pout.y();
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  LorentzSpinor<double> spt = sp.wave();
  complex<Energy> p0p3=pout.e()+pout.z();
  complex<Energy> p0m3=pout.e()-pout.z();
  s1 = UnitRemoval::InvE * 
    fact*( _left*mass*spt.s1()+_right*(p0m3*spt.s3()-p1m2*spt.s4()));
  s2 = UnitRemoval::InvE * 
    fact*( _left*mass*spt.s2()+_right*(p0p3*spt.s4()-p1p2*spt.s3()));
  s3 = UnitRemoval::InvE * 
    fact*(_right*mass*spt.s3()+ _left*(p0p3*spt.s1()+p1m2*spt.s2()));
  s4 = UnitRemoval::InvE * 
    fact*(_right*mass*spt.s4()+ _left*(p0m3*spt.s2()+p1p2*spt.s1()));
  return SpinorWaveFunction(pout,out,s1,s2,s3,s4);
}

// off-shell SpinorBar
SpinorBarWaveFunction FFSVertex::evaluate(Energy2 q2,int iopt,tcPDPtr out,
					  const SpinorBarWaveFunction & sbar,
					  const ScalarWaveFunction & sca,
					  complex<Energy> mass,
					  complex<Energy> width) {
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = sbar.momentum()+sca.momentum();
  // first calculate the couplings
  setCoupling(q2,out,sbar.particle(),sca.particle());
  Energy2 p2   = pout.m2();
  Complex fact = -norm()*sca.wave()*propagator(iopt,p2,out,mass,width);
  Complex ii(0.,1.);
  // momentum components
  if(mass.real() < ZERO) mass = out->mass();
  complex<Energy> p1p2 = pout.x()+ii*pout.y();
  complex<Energy> p1m2 = pout.x()-ii*pout.y();
  // complex numbers for the spinor
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  LorentzSpinorBar<double> sbart=sbar.wave();
  complex<Energy> p0p3=pout.e() +   pout.z();
  complex<Energy> p0m3=pout.e() -   pout.z();
  s1 = UnitRemoval::InvE * 
    fact*( mass*_left*sbart.s1()-_right*(p0p3*sbart.s3()+p1p2*sbart.s4()));
  s2 = UnitRemoval::InvE * 
    fact*( mass*_left*sbart.s2()-_right*(p1m2*sbart.s3()+p0m3*sbart.s4()));
  s3 = UnitRemoval::InvE * 
    fact*(mass*_right*sbart.s3()- _left*(p0m3*sbart.s1()-p1p2*sbart.s2()));
  s4 = UnitRemoval::InvE * 
    fact*(mass*_right*sbart.s4()+ _left*(p1m2*sbart.s1()-p0p3*sbart.s2()));
  return SpinorBarWaveFunction(pout,out,s1,s2,s3,s4);
}    
