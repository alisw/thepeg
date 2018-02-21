// -*- C++ -*-
//
// VVTVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVTVertex class.
//

#include "VVTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<VVTVertex> VVTVertex::initVVTVertex;
// Definition of the static class description member.

void VVTVertex::Init() {
    
  static ClassDocumentation<VVTVertex> documentation
    ("The VVTVertex class is the implementation of"
     "the vector-vector tensor vertices for helicity "
     "amplitude calculations. All such vertices should inherit"
     "from it.");
  
}

// function to evaluate the vertex
Complex VVTVertex::evaluate(Energy2 q2, const VectorWaveFunction & vec1,
			    const VectorWaveFunction & vec2, 
			    const TensorWaveFunction & ten,
			    Energy vmass) {
  // set the couplings
  setCoupling(q2,vec1.particle(),vec2.particle(),ten.particle());
  // mass of the vector
  if(vmass<ZERO) vmass = vec1.particle()->mass();
  // mass+k1.k2
  Energy2 mdot = vec1.momentum()*vec2.momentum();
  if(vmass!=ZERO) mdot += sqr(vmass);
  // dot product of wavefunctions and momenta
  Complex dotv1v2         = vec1.wave().dot(vec2.wave());
  complex<Energy> dotk1v2 = vec1.momentum()*vec2.wave();
  complex<Energy> dotk2v1 = vec1.wave()*vec2.momentum();
  // dot product of wavefunctions and momenta with the tensor
  LorentzPolarizationVector  tv1 = ten.wave().postDot(vec1.wave());
  LorentzPolarizationVector  tv2 = ten.wave().postDot(vec2.wave());
  LorentzPolarizationVectorE tk1 = ten.wave().postDot(vec1.momentum());
  LorentzPolarizationVectorE tk2 = ten.wave().postDot(vec2.momentum());
  Complex tenv1v2 = 
    vec1.wave    ().dot(tv2) + vec2.wave    ().dot(tv1);
  complex<Energy> tenk1v2 =
    vec1.momentum().dot(tv2) + vec2.wave    ().dot(tk1);
  complex<Energy> tenk2v1 =
    vec2.momentum().dot(tv1) + vec1.wave    ().dot(tk2);
  complex<Energy2> tenk1k2 =
    vec2.momentum().dot(tk1) + vec1.momentum().dot(tk2);
  // trace of the tensor
  Complex trace = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // evaluate the vertex
  return -0.5*Complex(0.,1.)*norm()*UnitRemoval::InvE2 *
    (trace*(dotk1v2*dotk2v1-dotv1v2*mdot)
     +mdot*tenv1v2-dotk2v1*tenk1v2
     -dotk1v2*tenk2v1+dotv1v2*tenk1k2);
}

// evaluate an off-shell vector
VectorWaveFunction VVTVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
				       const VectorWaveFunction & vec,
				       const TensorWaveFunction & ten,
				       complex<Energy> mass,
				       complex<Energy> width) {
  // evaluate the couplings
  setCoupling(q2,vec.particle(),out,ten.particle());
  // outgoing momentum
  Lorentz5Momentum pout = ten.momentum()+vec.momentum();
  // normalisation factor
  if(mass.real() < ZERO) mass   = out->mass();
  complex<Energy2> mass2 = sqr(mass);
  Energy2 p2    = pout.m2();
  Complex fact  = -0.5*norm()*propagator(iopt,p2,out,mass,width);
  // dot product of wavefunctions and momenta
  complex<Energy2> dotk1k2 = vec.momentum()*pout;
  complex<Energy> dotk2v1  = vec.wave()       *pout;
  // mass-k1.k2
  complex<Energy2> mdot = -dotk1k2;
  if(mass.real()!=ZERO) mdot += mass2;
  // components of the tensor
  Complex tentx = ten.tx()+ten.xt();
  Complex tenty = ten.ty()+ten.yt();
  Complex tentz = ten.tz()+ten.zt();
  Complex tenxy = ten.xy()+ten.yx();
  Complex tenxz = ten.xz()+ten.zx();
  Complex tenyz = ten.yz()+ten.zy();
  // dot product of momenta and polarization vectors with the tensor
  complex<Energy> tenk2v1 =
    2.*(+ten.tt()*vec.t()*pout.e() +ten.xx()*vec.x()*pout.x()
        +ten.yy()*vec.y()*pout.y()+ten.zz()*vec.z()*pout.z())
    -tentx*(vec.t()*pout.x()+vec.x()*pout.e())
    -tenty*(vec.t()*pout.y()+vec.y()*pout.e())
    -tentz*(vec.t()*pout.z()+vec.z()*pout.e())
    +tenxy*(vec.x()*pout.y()+vec.y()*pout.x())
    +tenxz*(vec.x()*pout.z()+vec.z()*pout.x())
    +tenyz*(vec.y()*pout.z()+vec.z()*pout.y());
  complex<Energy2> tenk1k2 =
    2.*(+ten.tt()*vec.e()*pout.e()  +ten.xx()*vec.px()*pout.x()
        +ten.yy()*vec.py()*pout.y()+ten.zz()*vec.pz()*pout.z())
    -tentx*(vec.e()*pout.x() +vec.px()*pout.e())
    -tenty*(vec.e()*pout.y() +vec.py()*pout.e())
    -tentz*(vec.e()*pout.z() +vec.pz()*pout.e())
    +tenxy*(vec.px()*pout.y()+vec.py()*pout.x())
    +tenxz*(vec.px()*pout.z()+vec.pz()*pout.x())
    +tenyz*(vec.py()*pout.z()+vec.pz()*pout.y());
  // trace of the tensor
  Complex trace = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // compute the vector
  Complex vec1[4];
  vec1[0] = UnitRemoval::InvE2*
    (mdot*(+   tentx*   vec.t()-2.*ten.xx()*vec.x()
	   -   tenxy*   vec.y()-    tenxz*  vec.z()-trace*vec.x())
     +(tenk2v1-trace*dotk2v1)*vec.px()-tenk1k2*vec.x()
     +dotk2v1*(+   tentx*   vec.e() -2.*ten.xx()*vec.px()
	       -   tenxy*   vec.py()-    tenxz*  vec.pz()));
  vec1[1] = UnitRemoval::InvE2*(
    mdot *
    (+    tenty*  vec.t()
     -    tenxy*  vec.x()
     -2.*ten.yy()*vec.y()
     -    tenyz*  vec.z()-trace*vec.y()
     )
    +(tenk2v1-trace*dotk2v1)*vec.py()
    -tenk1k2*vec.y()
    +dotk2v1*(+    tenty*  vec.e() 
	      -    tenxy*  vec.px()
	      -2.*ten.yy()*vec.py()
	      -    tenyz*  vec.pz()));

  
  vec1[2] = UnitRemoval::InvE2*
    (mdot*
     (+   tentz*   vec.t()-    tenxz*  vec.x()
      -   tenyz*   vec.y()-2.*ten.zz()*vec.z()-trace*vec.z())
     +(tenk2v1-trace*dotk2v1)*vec.pz()-tenk1k2*vec.z()
     +dotk2v1*(+   tentz*   vec.e() -    tenxz  *vec.px()
	       -   tenyz*   vec.py()-2.*ten.zz()*vec.pz()));
  vec1[3] = UnitRemoval::InvE2*
    (mdot*(+2.*ten.tt()*vec.t()-    tentx*  vec.x()
	   -   tenty*   vec.y()-    tentz*  vec.z()-trace*vec.t())
     +(tenk2v1-trace*dotk2v1)*vec.e() -tenk1k2*vec.t()
     +dotk2v1*(+2.*ten.tt()*vec.e() -    tentx*  vec.px()
	       -   tenty*   vec.py()-    tentz*  vec.pz()));
  // now add the piece for massive bosons
  if(mass.real()!=ZERO) {
    // DGRELL unit problem?
    Complex dot = tenk2v1 * UnitRemoval::InvE 
      -  dotk1k2 * trace  * UnitRemoval::InvE2;
    vec1[0] -= dot*pout.x() * UnitRemoval::InvE;
    vec1[1] -= dot*pout.y() * UnitRemoval::InvE;
    vec1[2] -= dot*pout.z() * UnitRemoval::InvE;
    vec1[3] -= dot*pout.e() * UnitRemoval::InvE;
  }
  // return the VectorWaveFunction
  for(int ix=0;ix<4;++ix){vec1[ix]=vec1[ix]*fact;}
  return VectorWaveFunction(pout,out,vec1[0],vec1[1],vec1[2],vec1[3]);
}

// off-shell tensor
TensorWaveFunction VVTVertex::evaluate(Energy2 q2, int iopt,tcPDPtr out,
				       const VectorWaveFunction & vec1,
				       const VectorWaveFunction & vec2,
				       Energy vmass,
				       complex<Energy> tmass,
				       complex<Energy> width) {
  // coupling
  setCoupling(q2,vec1.particle(),vec2.particle(),out);
  // momenta of the outgoing tensor
  // outgoing momentum
  Lorentz5Momentum pout= vec1.momentum()+vec2.momentum();
  // overall normalisation
  if(tmass.real() < ZERO) tmass   = out->mass();
  complex<Energy2> tmass2 = sqr(tmass);
  if(vmass<ZERO) vmass = vec1.particle()->mass();
  Energy2 vmass2 = sqr(vmass);
  Energy2 p2     = pout.m2();
  Complex fact   = 0.25*norm()*propagator(iopt,p2,out,tmass,width);
  // dot products we need to construct the tensor
  complex<Energy2> dotk1k2 = vec1.momentum()*vec2.momentum();
  complex<Energy> dotv1k2  = vec1.wave()*vec2.momentum();
  Complex dotv1v2          = vec1.wave().dot(vec2.wave());
  complex<Energy> dotk1v2  = vec1.momentum()*vec2.wave();
  complex<Energy2> dotkk1  = vec1.momentum()*pout;
  complex<Energy2> dotkk2  = vec2.momentum()*pout;
  complex<Energy> dotkv1   = vec1.wave()*pout;
  complex<Energy> dotkv2   = vec2.wave()*pout;
  // dot product ma^2+k1.k2
  complex<Energy2> mdot = vmass2+dotk1k2;
  // vectors to help construct the tensor
  Complex vecv1[4],vecv2[4];
  complex<Energy> veck1[4],veck2[4];
  complex<InvEnergy2> tmass2inv(ZERO);
  if(tmass.real()> ZERO) tmass2inv =  1./tmass2;
  vecv1[0]=vec1.x() -pout.x()*dotkv1*tmass2inv;
  vecv2[0]=vec2.x() -pout.x()*dotkv2*tmass2inv;
  veck1[0]=vec1.px()-pout.x()*dotkk1*tmass2inv;
  veck2[0]=vec2.px()-pout.x()*dotkk2*tmass2inv;
    
  vecv1[1]=vec1.y() -pout.y()*dotkv1*tmass2inv;
  vecv2[1]=vec2.y() -pout.y()*dotkv2*tmass2inv;
  veck1[1]=vec1.py()-pout.y()*dotkk1*tmass2inv;
  veck2[1]=vec2.py()-pout.y()*dotkk2*tmass2inv;
    
  vecv1[2]=vec1.z() -pout.z()*dotkv1*tmass2inv;
  vecv2[2]=vec2.z() -pout.z()*dotkv2*tmass2inv;
  veck1[2]=vec1.pz()-pout.z()*dotkk1*tmass2inv;
  veck2[2]=vec2.pz()-pout.z()*dotkk2*tmass2inv;
    
  vecv1[3]=vec1.t() -pout.e()*dotkv1*tmass2inv;
  vecv2[3]=vec2.t() -pout.e()*dotkv2*tmass2inv;
  veck1[3]=vec1.e() -pout.e()*dotkk1*tmass2inv;
  veck2[3]=vec2.e() -pout.e()*dotkk2*tmass2inv;
    
  // coefficient of g(nu,mu)-k^muk^nu/m^2
  Complex coeff1 = UnitRemoval::InvE2 * 
    (
     +4./3.*mdot*(-2.*dotv1v2+Complex((dotkv1*dotkv2+p2*dotv1v2)*tmass2inv))
     +4./3.*(dotv1k2*(dotk1v2-dotkk1*dotkv2*tmass2inv)
	     +dotk1v2*(dotv1k2-dotkk2*dotkv1*tmass2inv)
	     -dotv1v2*(dotk1k2-dotkk1*dotkk2*tmass2inv)
	     +(1.-p2*tmass2inv)*dotk1v2*dotv1k2)
     );
  // coefficient of g(nu,mu)
  Complex coeff2 = UnitRemoval::InvE2 * 
    (
     2.*mdot*(1.-p2*tmass2inv)*dotv1v2
     -2.*(1.-p2*tmass2inv)*dotk1v2*dotv1k2
     );
  // construct the tensor
  Complex ten[4][4];
  const Energy pout_tmp[4] = {pout.x(), pout.y(), pout.z(), pout.e()};
  for(int ix=0;ix<4;++ix) {
    for(int iy=0;iy<4;++iy) {
      complex<Energy2> temp;
      temp  = 2.*mdot*   (vecv1[ix]*vecv2[iy]+vecv1[iy]*vecv2[ix]);
      temp -= 2.*dotv1k2*(veck1[ix]*vecv2[iy]+veck1[iy]*vecv2[ix]);
      temp -= 2.*dotk1v2*(veck2[ix]*vecv1[iy]+veck2[iy]*vecv1[ix]);
      temp += 2.*dotv1v2*(veck1[ix]*veck2[iy]+veck1[iy]*veck2[ix]);
      
      ten[ix][iy] = UnitRemoval::InvE2 * temp
	-coeff1*tmass2inv*pout_tmp[ix]*pout_tmp[iy];
    }
  }
  // add the g(mu,nu) term
  ten[0][0] = ten[0][0]-(coeff1+coeff2);
  ten[1][1] = ten[1][1]-(coeff1+coeff2);
  ten[2][2] = ten[2][2]-(coeff1+coeff2);
  ten[3][3] = ten[3][3]+(coeff1+coeff2);
  // overall coefficent
  for(int ix=0;ix<4;++ix) {
    for(int iy=0;iy<4;++iy) {
      ten[ix][iy] = fact*ten[ix][iy];
    }
  }
  // return the wavefunction
  return TensorWaveFunction(pout,out,
			    ten[0][0],ten[0][1],ten[0][2],ten[0][3],
			    ten[1][0],ten[1][1],ten[1][2],ten[1][3],
			    ten[2][0],ten[2][1],ten[2][2],ten[2][3],
			    ten[3][0],ten[3][1],ten[3][2],ten[3][3]);
}
