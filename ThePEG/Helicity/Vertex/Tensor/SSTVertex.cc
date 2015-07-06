// -*- C++ -*-
//
// SSTVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSTVertex class.
//

#include "SSTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<SSTVertex> SSTVertex::initSSTVertex;
// Definition of the static class description member.

void SSTVertex::Init() {
  
  static ClassDocumentation<SSTVertex> documentation
("The SSTVertex class is the implementation of the "
 "helicity amplitude calculation for the scalar-scalar-tensor vertex"
 ", all such vertices should inherit from it");
  
}
// evaluate the vertex
Complex SSTVertex::evaluate(Energy2 q2, const ScalarWaveFunction & sca1,
    				const ScalarWaveFunction & sca2,
    				const TensorWaveFunction & ten) {
  // obtain the coupling
  setCoupling(q2,sca1.particle(),sca2.particle(),ten.particle());
  Complex ii(0.,1.);
  // evaluate the trace of the tensor
  Complex trace = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // dot product of the two momenta
  Energy2 dot = sca1.momentum()*sca2.momentum();
  Energy mass = sca1.particle()->mass();
  // second term
  complex<Energy2> second = 
    +2.*ten.tt()*sca1.e()*sca2.e()  +2.*ten.xx()*sca1.px()*sca2.px()
    +2.*ten.yy()*sca1.py()*sca2.py()+2.*ten.zz()*sca1.pz()*sca2.pz()
    -(ten.tx()+ten.xt())*(sca1.e()*sca2.px() +sca1.px()*sca2.e())
    -(ten.ty()+ten.yt())*(sca1.e()*sca2.py() +sca1.py()*sca2.e())
    -(ten.tz()+ten.zt())*(sca1.e()*sca2.pz() +sca1.pz()*sca2.e())
    +(ten.xy()+ten.yx())*(sca1.py()*sca2.px()+sca1.px()*sca2.py())
    +(ten.xz()+ten.zx())*(sca1.pz()*sca2.px()+sca1.px()*sca2.pz())
    +(ten.yz()+ten.zy())*(sca1.pz()*sca2.py()+sca1.py()*sca2.pz());
  // return the answer
  return -0.5*ii*norm()*UnitRemoval::InvE2*
    (trace*(mass*mass-dot)+second)*sca1.wave()*sca2.wave();
}
// off-shell tensor
TensorWaveFunction SSTVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
				       const ScalarWaveFunction & sca1,
				       const ScalarWaveFunction & sca2,
				       complex<Energy> mass, complex<Energy> width) {
  // obtain the coupling
  setCoupling(q2,sca1.particle(),sca2.particle(),out);
  // array for the tensor
  Complex ten[4][4];
  // calculate the outgoing momentum
  Lorentz5Momentum pout = sca1.momentum()+sca2.momentum();
  // prefactor
  if(mass.real() < ZERO) mass   = out->mass();
  complex<Energy2> mass2 = sqr(mass);
  Energy2 p2    = pout.m2();
  Complex fact  = 0.25*norm()*sca1.wave()*sca2.wave()*
    propagator(iopt,p2,out,mass,width);
  // dot products we need
  Energy2 dot12 = sca1.momentum()*sca2.momentum();
  Energy2 dot1  = sca1.momentum()*pout;
  Energy2 dot2  = pout*sca2.momentum();
  // the vectors that we need for the tensor
  LorentzPolarizationVectorE vec1,vec2;
  Complex a,b;
  Energy2 mphi2 = sqr(sca1.particle()->mass());
  // massive case
  if(mass.real()!=ZERO) {
    Complex norm1=dot1/mass2;
    Complex norm2=dot2/mass2;
    a = UnitRemoval::InvE2 * ((mphi2+dot12)*(Complex(2.*p2/mass2)-5.)
			      +4.*(dot12-dot1*dot2/mass2))/3.;
    b = -(-(mphi2+dot12)*(2.+p2/mass2)+4.*(dot12-dot1*(dot2/mass2)))/3./mass2;
    vec1 = LorentzPolarizationVectorE(sca1.momentum()) - 
      LorentzPolarizationVectorE(norm1 * pout);
    vec2 = LorentzPolarizationVectorE(sca2.momentum()) - 
      LorentzPolarizationVectorE(norm2 * pout);
  }
  // massless case
  else {
    a = UnitRemoval::InvE2 * (-5.*(mphi2+dot12)+4.*dot12)/3.;
    b = 0.;
    vec1 = sca1.momentum();
    vec2 = sca2.momentum();
  }
  // calculate the wavefunction
  complex<Energy> vec1_tmp[4] = {vec1.x(), vec1.y(), vec1.z(), vec1.t()};
  complex<Energy> vec2_tmp[4] = {vec2.x(), vec2.y(), vec2.z(), vec2.t()};
  complex<Energy> pout_tmp[4] = {pout.x(), pout.y(), pout.z(), pout.t()};
  for(int ix=0;ix<4;++ix) {
    for(int iy=0;iy<4;++iy) {
      complex<Energy2> temp = -2.*( vec1_tmp[ix]*vec2_tmp[iy] +
				    vec1_tmp[ix]*vec2_tmp[iy])
	-b*pout_tmp[ix]*pout_tmp[iy];
      ten[ix][iy]= UnitRemoval::InvE2 * temp;
    }
  }
  ten[3][3]=ten[3][3]-a;
  for(int ix=0;ix<3;++ix){ten[ix][ix]=ten[ix][ix]+a;}
  // prefactor
  for(int ix=0;ix<4;++ix){for(int iy=0;iy<4;++iy){ten[ix][iy]=fact*ten[ix][iy];}}
  // return the wavefunction
  return TensorWaveFunction(pout,out,
			    ten[0][0],ten[0][1],ten[0][2],ten[0][3],
			    ten[1][0],ten[1][1],ten[1][2],ten[1][3],
			    ten[2][0],ten[2][1],ten[2][2],ten[2][3],
			    ten[3][0],ten[3][1],ten[3][2],ten[3][3]);
}
// off-shell scalar
ScalarWaveFunction SSTVertex::evaluate(Energy2 q2,int iopt, tcPDPtr out,
				       const ScalarWaveFunction & sca,
				       const TensorWaveFunction & ten,
				       complex<Energy> mass, complex<Energy> width) {
  // obtain the coupling
  setCoupling(q2,sca.particle(),out,ten.particle());
  // calculate the outgoing momentum
  Lorentz5Momentum pout = sca.momentum()+ten.momentum();
  // prefactors
  if(mass.real() < ZERO) mass   = out->mass();
  complex<Energy2> mass2 = sqr(mass);
  Energy2 p2    = pout.m2();
  Complex fact  = 0.5*norm()*sca.wave()*propagator(iopt,p2,out,mass,width);
  // trace of the tensor
  Complex trace1 = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // dot product of the two momenta
  Energy2 dot = sca.momentum()*pout;
  // first term
  complex<Energy2> trace = trace1*(mass2-dot);
  // second term
  complex<Energy2> second = 
    +2.*ten.tt()*sca.e()*pout.e()  +2.*ten.xx()*sca.px()*pout.x()
    +2.*ten.yy()*sca.py()*pout.y()+2.*ten.zz()*sca.pz()*pout.z()
    -(ten.tx()+ten.xt())*( sca.e()*pout.x()+sca.px()*pout.e())
    -(ten.ty()+ten.yt())*( sca.e()*pout.y()+sca.py()*pout.e())
    -(ten.tz()+ten.zt())*( sca.e()*pout.z()+sca.pz()*pout.e())
    +(ten.xy()+ten.yx())*(sca.py()*pout.x()+sca.px()*pout.y())
    +(ten.xz()+ten.zx())*(sca.pz()*pout.x()+sca.px()*pout.z())
    +(ten.yz()+ten.zy())*(sca.py()*pout.z()+sca.pz()*pout.y());
  // put it all together
  second  = fact*(trace+second);
  Complex result = second * UnitRemoval::InvE2;
  return ScalarWaveFunction(pout,out,result);
}
