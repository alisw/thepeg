// -*- C++ -*-
//
// FFTVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFTVertex class.
//

#include "FFTVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// Definition of the static class description member
// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<FFTVertex,AbstractFFTVertex>
describeThePEGFFTVertex("ThePEG::FFTVertex", "libThePEG.so");
    
void FFTVertex::Init() {
  
  static ClassDocumentation<FFTVertex> documentation
    ("The FFTVertex class is the implementation of"
     "the fermion-antifermion tensor vertices for helicity "
     "amplitude calculations. All such vertices should inherit"
     "from it.");
  
}

// function to evaluate the vertex
Complex FFTVertex::evaluate(Energy2 q2,const SpinorWaveFunction & sp,
			    const SpinorBarWaveFunction & sbar,
			    const TensorWaveFunction & ten) {
  // set the couplings
  setCoupling(q2,sp.particle(),sbar.particle(),ten.particle());
  // vector current
  LorentzPolarizationVector as = sp.wave().vectorCurrent(sbar.wave());
  // momentum difference
  Lorentz5Momentum vdiff = sp.momentum()-sbar.momentum();
  // first term
  LorentzPolarizationVectorE test 
    = ten.wave().postDot(vdiff) + ten.wave().preDot(vdiff);
  complex<Energy> term1 = as.dot(test);
  // trace of polarization tensor
  Complex trace = ten.wave().trace();
  // dot products with polarization tensor
  // product of spinors
  Complex ffbar=  sp.wave().scalar(sbar.wave());
  // put everything together
  return -0.125*Complex(0.,1.)*norm()*UnitRemoval::InvE*
    ( term1 + 4.0*(sp.particle()->mass())*trace*ffbar );
}

// member function to evaluate an off-shell spinor
SpinorWaveFunction FFTVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
				       const SpinorWaveFunction & sp,
				       const TensorWaveFunction & ten,
				       complex<Energy> mass,
				       complex<Energy> width) {
  // momentum of the outgoing fermion
  Lorentz5Momentum pout = ten.momentum()+sp.momentum();
  // set the couplings
  setCoupling(q2,sp.particle(),out,ten.particle());
  Complex ii(0.,1.);
  // trace of the tensor
  Complex trace = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // mass of the fermion
  if(mass.real() < ZERO) mass = out->mass();
  // overall factor
  Energy2 p2 = pout.m2();
  Complex fact = 0.125*norm()*propagator(iopt,p2,out,mass,width);
  // compute the vector we need
  complex<Energy> dot[4];
  for(int ix=0;ix<4;++ix) {
    // evaluate the products we need
    dot[ix] =(ten(ix,3)+ten(3,ix))*(pout.e()+sp.e());
    dot[ix]-=(ten(ix,0)+ten(0,ix))*(pout.x()+sp.px());
    dot[ix]-=(ten(ix,1)+ten(1,ix))*(pout.y()+sp.py());
    dot[ix]-=(ten(ix,2)+ten(2,ix))*(pout.z()+sp.pz());
  }
  LorentzVector<complex<Energy> > vec(dot[0],dot[1],dot[2],dot[3]);
  vec -= 2.0 * trace * (pout + sp.momentum());
  // combinations of the vector
  complex<Energy> a1p2=vec.x()+ii*vec.y();
  complex<Energy> a1m2=vec.x()-ii*vec.y();
  LorentzSpinor<double>    spt  =sp  .wave();
  // now compute the first stage of the spinor wavefunction
  complex<Energy> a0p3=vec.t()+vec.z();
  complex<Energy> a0m3=vec.t()-vec.z();
  vec.setX(a0m3*spt.s3()-a1m2*spt.s4()); 
  vec.setY(a0p3*spt.s4()-a1p2*spt.s3());
  vec.setZ(a0p3*spt.s1()+a1m2*spt.s2());
  vec.setT(a0m3*spt.s2()+a1p2*spt.s1());
  if(mass.real()!=ZERO) {
    complex<Energy> dot = 4.*mass*trace;
    vec.setX(vec.x() + dot*spt.s1()); 
    vec.setY(vec.y() + dot*spt.s2());
    vec.setZ(vec.z() + dot*spt.s3());
    vec.setT(vec.t() + dot*spt.s4());
  }
  // combinations of the momentum
  complex<Energy> p1p2=pout.x()+ii*pout.y();
  complex<Energy> p1m2=pout.x()-ii*pout.y();
  // finally put everything together as the spinor
  Complex ferm[4];
  complex<Energy> p0p3=pout.e() +   pout.z();
  complex<Energy> p0m3=pout.e() -   pout.z();
  ferm[0] = Complex(UnitRemoval::InvE2 * fact*( p0m3*vec.z()-p1m2*vec.t()));
  ferm[1] = Complex(UnitRemoval::InvE2 * fact*(-p1p2*vec.z()+p0p3*vec.t()));
  ferm[2] = Complex(UnitRemoval::InvE2 * fact*( p0p3*vec.x()+p1m2*vec.y()));
  ferm[3] = Complex(UnitRemoval::InvE2 * fact*( p1p2*vec.x()+p0m3*vec.y()));
  if(mass.real()!=ZERO) {
    ferm[0] += Complex(UnitRemoval::InvE2 * fact*(mass*vec.x()));
    ferm[1] += Complex(UnitRemoval::InvE2 * fact*(mass*vec.y()));
    ferm[2] += Complex(UnitRemoval::InvE2 * fact*(mass*vec.z()));
    ferm[3] += Complex(UnitRemoval::InvE2 * fact*(mass*vec.t()));
  }
  // return the wavefunction
  return SpinorWaveFunction(pout,out,ferm[0],ferm[1],ferm[2],ferm[3]);
}


// member function to evaluate an off-shell spinor bar
SpinorBarWaveFunction FFTVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
					  const SpinorBarWaveFunction & sbar,
					  const TensorWaveFunction & ten,
					  complex<Energy> mass,
					  complex<Energy> width) {
  // momentum of the outgoing fermion
  Lorentz5Momentum pout = ten.momentum()+sbar.momentum();
  // set the couplings
  setCoupling(q2,out,sbar.particle(),ten.particle());
  Complex ii(0.,1.);
  // trace of the tensor
  Complex trace = ten.tt()-ten.xx()-ten.yy()-ten.zz();
  // mass of the fermion
  if(mass.real() < ZERO) mass = out->mass();
      // overall factor
  Energy2 p2 = pout.m2();
  Complex fact=0.125*norm()*propagator(iopt,p2,out,mass,width);
  // vector
  complex<Energy> dot[4];
  for(int ix=0;ix<4;++ix) {
    // evaluate the products we need
    dot[ix] =-(ten(ix,3)+ten(3,ix))*(pout.e()+sbar.e());
    dot[ix]+= (ten(ix,0)+ten(0,ix))*(pout.x()+sbar.px());
    dot[ix]+= (ten(ix,1)+ten(1,ix))*(pout.y()+sbar.py());
    dot[ix]+= (ten(ix,2)+ten(2,ix))*(pout.z()+sbar.pz());
  }
  LorentzVector<complex<Energy> > vec(dot[0],dot[1],dot[2],dot[3]);
  vec += 2.*trace*(pout+sbar.momentum());
  // combinations of the vector
  complex<Energy> a1p2=vec.x()+ii*vec.y();
  complex<Energy> a1m2=vec.x()-ii*vec.y();
  LorentzSpinorBar<double> sbart=sbar.wave();
  // now compute the first stage of the spinorbar wavefunction
  complex<Energy> a0p3=vec.t()+vec.z();
  complex<Energy> a0m3=vec.t()-vec.z();
  vec.setX(a0p3*sbart.s3()+a1p2*sbart.s4()); 
  vec.setY(a0m3*sbart.s4()+a1m2*sbart.s3());
  vec.setZ(a0m3*sbart.s1()-a1p2*sbart.s2());
  vec.setT(a0p3*sbart.s2()-a1m2*sbart.s1());
  if(mass.real()!=ZERO) {
    complex<Energy> dot = 4.*mass*trace;
    vec.setX(vec.x() + dot*sbart.s1()); 
    vec.setY(vec.y() + dot*sbart.s2());
    vec.setZ(vec.z() + dot*sbart.s3());
    vec.setT(vec.t() + dot*sbart.s4());
  }
  // combinations of the momentum
  complex<Energy>  p1p2=pout.x()+ii*pout.y();
  complex<Energy>  p1m2=pout.x()-ii*pout.y();
  // finally put everything together as the spinor
  Complex ferm[4];
  complex<Energy> p0p3=pout.e() +   pout.z();
  complex<Energy> p0m3=pout.e() -   pout.z();
  ferm[0] = Complex(UnitRemoval::InvE2 * fact*(-p0p3*vec.z()-p1p2*vec.t()));
  ferm[1] = Complex(UnitRemoval::InvE2 * fact*(-p1m2*vec.z()-p0m3*vec.t()));
  ferm[2] = Complex(UnitRemoval::InvE2 * fact*(-p0m3*vec.x()+p1p2*vec.y()));
  ferm[3] = Complex(UnitRemoval::InvE2 * fact*( p1m2*vec.x()-p0p3*vec.y()));
  if(mass.real()!=ZERO) {
    ferm[0] += Complex(UnitRemoval::InvE2 * fact*mass*vec.x());
    ferm[1] += Complex(UnitRemoval::InvE2 * fact*mass*vec.y());
    ferm[2] += Complex(UnitRemoval::InvE2 * fact*mass*vec.z());
    ferm[3] += Complex(UnitRemoval::InvE2 * fact*mass*vec.t());
  }
  // return the wavefunction
  return SpinorBarWaveFunction(pout,out,ferm[0],ferm[1],ferm[2],ferm[3]);
}

// member function to evaluate an off-shell tensor
TensorWaveFunction FFTVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
				       const SpinorWaveFunction & sp,
				       const SpinorBarWaveFunction & sbar,
				       complex<Energy> mass,
				       complex<Energy> width) {
  // calculating the couplings
  setCoupling(q2,sp.particle(),sbar.particle(),out);
  Complex ii(0.,1.);
  // momentum of the outgoing tensor
  Lorentz5Momentum pout = sp.momentum()+sbar.momentum();
  // calculate the prefactor
  Energy2 p2   = pout.m2();
  Complex fact = 0.0625*norm()*propagator(iopt,p2,out,mass,width);
  if(mass.real() < ZERO) mass  = out->mass();
  complex<Energy2> mass2 = sqr(mass);
  // spinor vector
  Complex aspin[4];
  LorentzSpinorBar<double> sbart=sbar.wave();
  LorentzSpinor<double>    spt  =sp  .wave();
  aspin[3] = sbart.s1()*spt.s3()+sbart.s2()*spt.s4()
            +sbart.s3()*spt.s1()+sbart.s4()*spt.s2();
  // spatial components are the same in both conventions
  aspin[0] =     +sbart.s1()*spt.s4()+sbart.s2()*spt.s3()
                 -sbart.s3()*spt.s2()-sbart.s4()*spt.s1();
  aspin[1] = ii*(-sbart.s1()*spt.s4()+sbart.s2()*spt.s3()
		 +sbart.s3()*spt.s2()-sbart.s4()*spt.s1());
  aspin[2] =     +sbart.s1()*spt.s3()-sbart.s2()*spt.s4()
                 -sbart.s3()*spt.s1()+sbart.s4()*spt.s2();
  // mass dependent term
  Complex ffbar;
  if(sp.particle()->mass()!=ZERO) {
    ffbar = Complex(UnitRemoval::InvE * (sp.particle()->mass())*
		    (sp.s1()*sbar.s1()+sp.s2()*sbar.s2()+sp.s3()*sbar.s3()+sp.s4()*sbar.s4()));
  }
  else {
    ffbar = 0.;
  }
  // dot products for the calculation
  complex<Energy> dotka = 
    +aspin[3]*pout.e()-aspin[0]*pout.x()
    -aspin[1]*pout.y()-aspin[2]*pout.z();
  complex<Energy> dot12a = 
    +aspin[3]*(sp.e() -sbar.e() )-aspin[0]*(sp.px()-sbar.px())
    -aspin[1]*(sp.py()-sbar.py())-aspin[2]*(sp.pz()-sbar.pz());
  complex<Energy2> diff=sp.m2()-sbar.m2();
  complex<InvEnergy> dotkam=dotka/mass2;
  Complex diffm =diff/mass2;
  Complex p2m = p2/mass2;
  // construct the vectors for the first two terms
  Complex veca[4],vecb[4];

  veca[0] = aspin[0]-UnitRemoval::InvE2*dotka*pout.x();
  vecb[0] = Complex(UnitRemoval::InvE*(sp.px()-sbar.px()-diffm*pout.x()));
  
  veca[1] = aspin[1]-UnitRemoval::InvE2*dotka*pout.y();
  vecb[1] = Complex(UnitRemoval::InvE*(sp.py()-sbar.py()-diffm*pout.y()));
  
  veca[2] = aspin[2]-UnitRemoval::InvE2*dotka*pout.z();
  vecb[2] = Complex(UnitRemoval::InvE*(sp.pz()-sbar.pz()-diffm*pout.z()));
  
  veca[3] = aspin[3]-UnitRemoval::InvE2*dotka*pout.e();
  vecb[3] = Complex(UnitRemoval::InvE*(sp.e()-sbar.e()-diffm*pout.e()));
  
  // coefficients fr hte second two terms
  Complex temp = UnitRemoval::InvE*(p2m*dot12a-dotkam*diff);
  Complex coeff1 = -4./3.*(2.*ffbar*(1.-p2m) + temp);
  temp = Complex(UnitRemoval::InvE*(-3.*dot12a+2.*p2m*dot12a+diffm*dotka));
  Complex coeff2 = -4./3./mass2*( 4.*ffbar*(1.-p2m) + temp)*UnitRemoval::E2;
  // construct the tensor
  Complex ten[4][4];
  
  const complex<Energy> pout_tmp[4] 
    = {pout.x(), pout.y(), pout.z(), pout.e()};

  for(int ix=0;ix<4;++ix) {
    for(int iy=0;iy<4;++iy) {
      Complex temp = coeff2*pout_tmp[ix]*pout_tmp[iy]*UnitRemoval::InvE2;
      ten[ix][iy] = 2.*(veca[ix]*vecb[iy]+veca[iy]*vecb[ix]) + temp;
    }
  }
  
  ten[0][0]=ten[0][0]-coeff1;
  ten[1][1]=ten[1][1]-coeff1;
  ten[2][2]=ten[2][2]-coeff1;
  ten[3][3]=ten[3][3]+coeff1;
  // multiply by final prefactor
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
