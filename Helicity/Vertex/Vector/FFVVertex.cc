// -*- C++ -*-
//
// FFVVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFVVertex class.
//

#include "FFVVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// Definition of the static class description member
// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<FFVVertex,AbstractFFVVertex>
describeThePEGFFVVertex("ThePEG::FFVVertex", "libThePEG.so");
    
void FFVVertex::Init() {
      
  static ClassDocumentation<FFVVertex> documentation
    ("The FFVVertex class implements the helicity amplitude"
     "calculations for a fermion-fantifermion gauge boson vertex. Any   "
     "implementation of such a vertex should inherit from in and implement"
     " the virtual setCoupling member to calculate the coupling");
}

// evalulate the full vertex
Complex FFVVertex::evaluate(Energy2 q2,
			    const SpinorWaveFunction & sp, 
			    const SpinorBarWaveFunction & sbar,
			    const VectorWaveFunction & vec) {
  // first calculate the couplings
  if(kinematics()) calculateKinematics(sp.momentum(),sbar.momentum(),vec.momentum());
  setCoupling(q2,sp.particle(),sbar.particle(),vec.particle());
  Complex ii(0.,1.);
  // useful combinations of the polarization vector components
  Complex e0p3=vec.t()+vec.z();
  Complex e0m3=vec.t()-vec.z();
  Complex e1p2=vec.x()+ii*vec.y();
  Complex e1m2=vec.x()-ii*vec.y();
  Complex vertex(0.);
  if(_left!=0.) {
    vertex += _left*(+sbar.s3()*(sp.s1()*e0p3+sp.s2()*e1m2)
		     +sbar.s4()*(sp.s1()*e1p2+sp.s2()*e0m3)); 
  }
  // then the right piece (often not needed eg W vertex)
  if(_right!=0.) {
    vertex += _right*(+sbar.s1()*(sp.s3()*e0m3-sp.s4()*e1m2)
		      -sbar.s2()*(sp.s3()*e1p2-sp.s4()*e0p3));
  }
  vertex*=ii;
  return vertex*norm();
}

// evaluate an off-shell spinor
SpinorWaveFunction FFVVertex::evaluate(Energy2 q2, int iopt,tcPDPtr  out,
				       const SpinorWaveFunction & sp,
				       const VectorWaveFunction &vec,
				       complex<Energy> mass,
				       complex<Energy> width) {
  // extract the pointers to the particle data objects
  tcPDPtr  Psp=sp.particle();
  tcPDPtr  Pvec=vec.particle();
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = sp.momentum()+vec.momentum(); 
  // first calculate the couplings
  if(kinematics()) calculateKinematics(sp.momentum(),pout,vec.momentum());
  setCoupling(q2,Psp,out,Pvec);
  Complex ii(0.,1.);  // now evaluate the contribution
  // polarization components
  Complex e0p3 = vec.t() +  vec.z();
  Complex e0m3 = vec.t() -  vec.z();
  Complex e1p2 = vec.x()+ii*vec.y();
  Complex e1m2 = vec.x()-ii*vec.y();
  // overall factor
  Energy2 p2 = pout.m2();
  Complex fact=-normPropagator(iopt,p2,out,mass,width);
  // momentum components
  if(mass.real() < ZERO) mass  = iopt==5 ? ZERO : out->mass();
  complex<Energy> p1p2 = pout.x()+ii*pout.y();
  complex<Energy> p1m2 = pout.x()-ii*pout.y();
  // complex nos for for the spinor
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  LorentzSpinor<double>    spt  =sp  .wave();
  complex<Energy> p0p3=pout.e() +   pout.z();
  complex<Energy> p0m3=pout.e() -   pout.z();
  // left piece
  if(_left!=0.) {
    Complex a3=_left*fact*( spt.s1()*e0p3+spt.s2()*e1m2);
    Complex a4=_left*fact*( spt.s1()*e1p2+spt.s2()*e0m3);
    s1 += Complex(UnitRemoval::InvE * (p0m3*a3-p1m2*a4));
    s2 += Complex(UnitRemoval::InvE * (-p1p2*a3+p0p3*a4));
    s3 += Complex(UnitRemoval::InvE * a3*mass);
    s4 += Complex(UnitRemoval::InvE * a4*mass);
  }
  // right piece
  if(_right!=0.) {
    Complex a1=_right*fact*( spt.s3()*e0m3-spt.s4()*e1m2);
    Complex a2=_right*fact*(-spt.s3()*e1p2+spt.s4()*e0p3);
    s1 += Complex(UnitRemoval::InvE * a1*mass);
    s2 += Complex(UnitRemoval::InvE * a2*mass);
    s3 += Complex(UnitRemoval::InvE * (p0p3*a1+p1m2*a2));
    s4 += Complex(UnitRemoval::InvE * (p1p2*a1+p0m3*a2));
  }
  // return the wavefunction
  return SpinorWaveFunction(pout,out,s1,s2,s3,s4);
}
  
// evaluate an off-shell SpinorBar
SpinorBarWaveFunction FFVVertex::evaluate(Energy2 q2,int iopt,tcPDPtr  out,
					  const SpinorBarWaveFunction & sbar,
					  const VectorWaveFunction & vec,
					  complex<Energy> mass,
					  complex<Energy> width) {
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = sbar.momentum()+vec.momentum();
  // first calculate the couplings
  if(kinematics()) calculateKinematics(pout,sbar.momentum(),vec.momentum());
  setCoupling(q2,out,sbar.particle(),vec.particle());
  Complex ii(0.,1.);
  // now evaluate the contribution
  // polarization components
  Complex e0p3=vec.t() +  vec.z();
  Complex e0m3=vec.t() -  vec.z();
  Complex e1p2=vec.x()+ii*vec.y();
  Complex e1m2=vec.x()-ii*vec.y();
  // overall factor
  Energy2 p2 = pout.m2();
  Complex fact=-normPropagator(iopt,p2,out,mass,width);
  // momentum components
  if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
  complex<Energy> p1p2=pout.x()+ii*pout.y();
  complex<Energy> p1m2=pout.x()-ii*pout.y();
  // complex numbers for the spinor
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  LorentzSpinorBar<double> sbart=sbar.wave();
  Energy p0p3=pout.e() +   pout.z();
  Energy p0m3=pout.e() -   pout.z();
  // left piece
  if(_left!=0.) {
    Complex a1 =  _left*fact*( sbart.s3()*e0p3+sbart.s4()*e1p2);
    Complex a2 =  _left*fact*( sbart.s3()*e1m2+sbart.s4()*e0m3);
    s1 += Complex(UnitRemoval::InvE*a1*mass);
    s2 += Complex(UnitRemoval::InvE*a2*mass);
    s3 += Complex(UnitRemoval::InvE*(-p0m3*a1+p1p2*a2));
    s4 += Complex(UnitRemoval::InvE*(+p1m2*a1-p0p3*a2));
  }
  // right piece
  if(_right!=0.) {
    Complex a3 = _right*fact*( sbart.s1()*e0m3-sbart.s2()*e1p2);
    Complex a4 = _right*fact*(-sbart.s1()*e1m2+sbart.s2()*e0p3);
    s1 += Complex(UnitRemoval::InvE*(-p0p3*a3-p1p2*a4));
    s2 += Complex(UnitRemoval::InvE*(-p1m2*a3-p0m3*a4));
    s3 += Complex(UnitRemoval::InvE*a3*mass);
    s4 += Complex(UnitRemoval::InvE*a4*mass);
  }
  return SpinorBarWaveFunction(pout,out,s1,s2,s3,s4);
}

// off-shell vector
VectorWaveFunction FFVVertex::evaluate(Energy2 q2,int iopt,tcPDPtr  out,
				       const SpinorWaveFunction & sp,
				       const SpinorBarWaveFunction & sbar,
				       complex<Energy> mass,
				       complex<Energy> width) {
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = sbar.momentum()+sp.momentum();
  // first calculate the couplings
  if(kinematics()) calculateKinematics(sp.momentum(),sbar.momentum(),pout);
  setCoupling(q2,sp.particle(),sbar.particle(),out);
  Complex ii(0.,1.);
  // overall factor
  Energy2 p2 = pout.m2();
  Complex fact = normPropagator(iopt,p2,out,mass,width);
  // momentum components
  if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
  complex<Energy2> mass2 = sqr(mass);
  // the vector for the fermion-antifermion
  Complex vec[4];
  // left coupling
  if(_left!=0.) {
    vec[0] =   -_left*(sbar.s3()*sp.s2()+sbar.s4()*sp.s1());
    vec[1] = ii*_left*(sbar.s3()*sp.s2()-sbar.s4()*sp.s1());
    vec[2] =   -_left*(sbar.s3()*sp.s1()-sbar.s4()*sp.s2());
    vec[3] =    _left*(sbar.s3()*sp.s1()+sbar.s4()*sp.s2());
  }
  // right coupling
  if(_right!=0.) {
    vec[0] +=    +_right*(sbar.s1()*sp.s4()+sbar.s2()*sp.s3());
    vec[1] += -ii*_right*(sbar.s1()*sp.s4()-sbar.s2()*sp.s3());
    vec[2] +=    +_right*(sbar.s1()*sp.s3()-sbar.s2()*sp.s4());
    vec[3] +=    +_right*(sbar.s1()*sp.s3()+sbar.s2()*sp.s4());
  }
  // massless boson
  if(mass.real()==ZERO) {
    for(int ix=0;ix<4;++ix){vec[ix]*=fact;}
  }
  // massive boson
  else {
    complex<InvEnergy> dot = ( pout.e() *vec[3]
			       -pout.x()*vec[0]
			       -pout.y()*vec[1]
			       -pout.z()*vec[2])/mass2;
    vec[0]=fact*(vec[0]-dot*pout.x());
    vec[1]=fact*(vec[1]-dot*pout.y());
    vec[2]=fact*(vec[2]-dot*pout.z());
    vec[3]=fact*(vec[3]-dot*pout.e());
  }
  return VectorWaveFunction(pout,out,vec[0],vec[1],vec[2],vec[3]);
}

SpinorWaveFunction FFVVertex::evaluateSmall(Energy2 q2,int iopt, tcPDPtr out,
					    const SpinorWaveFunction & sp,
					    const VectorWaveFunction & vec,
					    unsigned int fhel, unsigned int vhel,
					    double ctheta, double phi, double stheta,
					    bool includeEikonal,
					    SmallAngleDirection direction,
					    Energy mass, Energy) {
  assert(fhel <= 1);
  assert( vhel == 0 || vhel == 2 );
  SpinorWaveFunction output;
  // first calculate the couplings
  setCoupling(q2,sp.particle(),out,vec.particle());
  Complex ii(0.,1.);
  if(mass < ZERO) mass = iopt==5 ? ZERO : out->mass();
  Lorentz5Momentum pout = sp.momentum()+vec.momentum();
  assert(sp.direction()!=intermediate);
  // helicity of the boson
  double lam = double(vhel)-1.;
  // energy of the boson
  Energy Eg = abs(vec.momentum().e());
  // energy of the fermion
  Energy Ef = abs(sp .momentum().e());
  // energy fraction of photon
  double x = Eg/Ef;
  // velocity of the fermon
  double beta = sqrt(1.-sqr(mass/Ef));
  // dimensionless versions of the variables
  double dm  = mass*UnitRemoval::InvE;
  double dEf =   Ef*UnitRemoval::InvE;
  double rE = sqrt(.5*dEf);
  // calculation of propagator accurate as beta->1 and theta -> 0
  Energy2 dot = 2.*Ef*Eg*(sqr(mass/Ef)/(1.+beta)*ctheta 
			  + sqr(stheta)/(1.+ctheta) );
  Complex fact=  norm()*(0.5*left()+0.5*right())*UnitRemoval::E2/dot;
  // phase factor
  Complex ephig = cos(phi )+ii*sin(phi );
  // calculation of the spinor
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  // u-type spinor
  if(sp.wave().Type()==SpinorType::u) {
    // fermion along +z
    if(direction==PostiveZDirection) {
      if(fhel==0) {
	s1 = +x*rE*dEf*sqrt(1.+beta)/ephig*(1.+lam)*sqr(stheta)/(1.+ctheta);
	s2 =-rE*dEf*sqrt(1.+beta)*stheta*
	  (+x*(1.+lam)-(includeEikonal ? 2.*beta*lam : 0. ));
	s3 = +x*rE*dm/ephig*(lam-1.)*(1.+ctheta)/sqrt(1.+beta);
	s4 = +rE*dm/sqrt(1.+beta)*stheta*
	   (-x*(1.-lam)+(includeEikonal ? 2.*beta*lam : 0. ));
      }
      else if(fhel==1) {     
	s1 = +rE*dm/sqrt(1.+beta)*stheta*
	  (+x*(1.+lam)+(includeEikonal ? 2.*beta*lam : 0. ));
	s2 = -x*rE*dm*ephig*(lam+1.)*(1.+ctheta)/sqrt(1.+beta);
	s3 = -rE*dEf*sqrt(1.+beta)*stheta*
	  (-x*(1.-lam)-(includeEikonal ? 2.*beta*lam : 0. ));
	s4 = -x*rE*dEf*sqrt(1.+beta)*ephig*(lam-1.)*sqr(stheta)/(1.+ctheta);
      }
    }
    // fermion along -z
    else {
      if(fhel==0) {
	s1 = -rE*dEf*sqrt(1.+beta)*stheta*
	  (+x*(1.+lam)-(includeEikonal ? 2.*beta*lam : 0. ));
	s2 = -x*rE*dEf*sqrt(1.+beta)/ephig*(1.+lam)*sqr(stheta)/(1.+ctheta);
	s3 = rE*dm/sqrt(1.+beta)*stheta*
	  (-x*(1.-lam)+(includeEikonal ? 2.*beta*lam : 0. ));
	s4 = -x*rE*dm/ephig*(-1.+lam)*(1.+ctheta)/sqrt(1.+beta);
      }
      else if(fhel==1) { 
	s1 =-x*rE*dm*ephig*(1.+lam)*(1.+ctheta)/sqrt(1.+beta);
	s2 =-rE*dm/sqrt(1.+beta)*stheta*
	  ( x*(1.+lam) + (includeEikonal ? 2.*beta*lam : 0. ));
	s3 = x*rE*dEf*sqrt(1.+beta)*ephig*(1.-lam)*sqr(stheta)/(1.+ctheta);
	s4 =-rE*dEf*sqrt(1.+beta)*stheta*
	  ( x*(1.-lam) + (includeEikonal ? 2.*beta*lam : 0. ));
      }
    }
  }
  // v-type spinor
  else if(sp.wave().Type()==SpinorType::v) {
    // fermion along +z
    if(direction==PostiveZDirection) {
      if(fhel==0) {
	s1 = rE*dm/sqrt(1.+beta)*stheta*
	  (-x*(1.+lam) + ( includeEikonal ? 2.*beta*lam : 0. ));
	s2 = x*rE*dm*ephig*(lam+1.)*(1.+ctheta)/sqrt(1.+beta);
	s3 = rE*dEf*sqrt(1.+beta)*stheta*
	  ( x*(1.-lam) - ( includeEikonal ? 2.*beta*lam : 0. ));
	s4 = x*rE*dEf*sqrt(1.+beta)*ephig*(1.-lam)*sqr(stheta)/(1.+ctheta);
      }
      else if(fhel==1) {
	s1 = x*rE*dEf*sqrt(1.+beta)/ephig*(1.+lam)*sqr(stheta)/(1.+ctheta);
	s2 =-rE*dEf*sqrt(1.+beta)*stheta*
	  ( x*(1.+lam) + (includeEikonal ? 2.*beta*lam : 0.));
	s3 = x*rE*dm/ephig*(1.-lam)*(1.+ctheta)/sqrt(1.+beta);
	s4 = rE*dm/sqrt(1.+beta)*stheta*
	  ( x*(1.-lam) + (includeEikonal ? 2.*beta*lam : 0. ));
      }
    }
    // fermion aling -z
    else {
      if(fhel==0) {
	s1 = x*rE*dm*ephig*(1.+lam)*(1.+ctheta)/sqrt(1.+beta);
	s2 = rE*dm/sqrt(1.+beta)*stheta*
	  ( x*(1.+lam) - ( includeEikonal ? 2.*beta*lam : 0. ));
	s3 = x*rE*dEf*sqrt(1.+beta)*ephig*(1.-lam)*sqr(stheta) / (1.+ctheta);
	s4 =-rE*dEf*sqrt(1.+beta)*stheta*
	  ( x*(1.-lam) - ( includeEikonal ? 2.*beta*lam : 0. ));

      }
      else if(fhel==1) {
	s1 =-rE*dEf*sqrt(1.+beta)*stheta*
	  ( x*(1.+lam) + ( includeEikonal ? 2.*beta*lam : 0. ));
	s2 =-x*rE*dEf*sqrt(1.+beta)/ephig*(1.+lam)*sqr(stheta)/(1.+ctheta);
	s3 = rE*dm/sqrt(1.+beta)*stheta*
	  ( x*(1.-lam) + ( includeEikonal ? 2.*beta*lam : 0. ));
	s4 =-x*rE*dm/ephig*(1.-lam)*(1.+ctheta)/sqrt(1.+beta);
      }
    }
  }
  s1 *= -fact;
  s2 *= -fact;
  s3 *= -fact;
  s4 *= -fact;  
  return SpinorWaveFunction(pout,out,s1,s2,s3,s4);
}


SpinorBarWaveFunction FFVVertex::evaluateSmall(Energy2 q2,int iopt, tcPDPtr out,
					       const SpinorBarWaveFunction & sbar,
					       const VectorWaveFunction & vec,
					       unsigned int fhel, unsigned int vhel,
					       double ctheta, double phi, double stheta,
					       bool includeEikonal,
					       SmallAngleDirection direction,
					       Energy mass, Energy) {
  assert(fhel <= 1);
  assert( vhel == 0 || vhel == 2 );
  SpinorBarWaveFunction output;
  // first calculate the couplings
  setCoupling(q2,out,sbar.particle(),vec.particle());
  Complex ii(0.,1.);
  if(mass < ZERO) mass  = iopt==5 ? ZERO : out->mass();
  Lorentz5Momentum pout = sbar.momentum()+vec.momentum();
  assert(sbar.direction()!=intermediate);
  // helicity of the boson
  double lam = double(vhel)-1.;
  // energies and velocities
  Energy Ef = abs(sbar.momentum().e());
  Energy Eg = abs(vec .momentum().e());
  double x = Eg/Ef;
  double beta = sqrt(1.-sqr(mass/Ef));
  // calculation of propagator accurate as beta->1 and theta -> 0
  Energy2 dot = 2.*Ef*Eg*(sqr(mass/Ef)/(1.+beta)*ctheta 
			  + sqr(stheta)/(1.+ctheta) );
  Complex fact=  norm()*(0.5*left()+0.5*right())*UnitRemoval::E2/dot;
  // calculation of the spinor
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  Complex ephig = cos(phi )+ii*sin(phi );
  double dm  = mass*UnitRemoval::InvE;
  double dEf =   Ef*UnitRemoval::InvE;
  double rE = sqrt(.5*dEf);
  // u-type spinor
  if(sbar.wave().Type()==SpinorType::u) {
    // fermion along +z
    if(direction==PostiveZDirection) {
      if(fhel==0) {    
	s1 = x*rE*dm*ephig*(1.+lam)*(1.+ctheta)/sqrt(1.+beta);
	s2 = rE*dm/sqrt(1.+beta)*stheta*
	  (+x*(1.+lam) - (includeEikonal ? 2.*beta*lam : 0. ));
	s3 = -x*rE*dEf*sqrt(1.+beta)*ephig*(1.-lam)*sqr(stheta)/(1.+ctheta);
	s4 = rE*dEf*sqrt(1.+beta)*stheta*
	  (+x*(1.-lam) - (includeEikonal ? 2.*beta*lam : 0. ));
      }
      else if(fhel==1) {
	s1 = -rE*dEf*sqrt(1.+beta)*stheta*
	  ( x*(1.+lam) + (includeEikonal ? 2.*beta*lam : 0. ));
	s2 = -x*rE*dEf*sqrt(1.+beta)/ephig*(lam+1.)*sqr(stheta)/(1.+ctheta);
	s3 =-rE*dm/sqrt(1.+beta)*stheta*
	  ( x*(1.-lam) + (includeEikonal ? 2.*beta*lam : 0. ));
	s4 = -x*rE*dm/ephig*(lam-1.)*(1.+ctheta)/sqrt(1.+beta);


	s4 = rE*dm*(1.+ctheta)/ephig*x*(1.-lam)/sqrt(1.+beta);
      }
    }
    // fermion aling -z
    else {
      if(fhel==0) {
	s1 =  rE*dm/sqrt(1.+beta)*stheta*
	  (+x*(1.+lam) - (includeEikonal ? 2.*beta*lam : 0. ));
	s2 = -x*rE*dm*ephig*(1.+lam)*(1.+ctheta)/sqrt(1.+beta);
	s3 = rE*dEf*sqrt(1.+beta)*stheta*
	  ( x*(1.-lam) - (includeEikonal ? 2.*beta*lam : 0. ));
	s4 = -x*rE*dEf*sqrt(1.+beta)*ephig*(lam-1.)*sqr(stheta)/(1.+ctheta);
      }
      else if(fhel==1) {
	s1 = -x*rE*dEf*sqrt(1.+beta)/ephig*(1.+lam)*sqr(stheta)/(1.+ctheta);
	s2 = rE*dEf*sqrt(1.+beta)*stheta*
	  (+x*(1.+lam) + (includeEikonal ? 2.*beta*lam : 0. ));
	s3 =-x*rE*dm/ephig*(lam-1.)*(1.+ctheta)/sqrt(1.+beta);
	s4 = rE*dm/sqrt(1.+beta)*stheta*
	  ( x*(1.-lam) + (includeEikonal ? 2.*beta*lam : 0. ));
      }
    }
  }
  // v-type spinor
  else if(sbar.wave().Type()==SpinorType::v) {
    // anti fermion along +z
    if(direction==PostiveZDirection) {
      if(fhel==0) {    
	s1 = -rE*dEf*sqrt(1.+beta)*stheta*
	  ( x*(1.+lam) - (includeEikonal ?  2.*beta*lam : 0. ));
	s2 = -x*rE*dEf*sqrt(1.+beta)/ephig*(1.+lam)*sqr(stheta)/(1.+ctheta);
	s3 = rE*dm/sqrt(1.+beta)*stheta*
	  ( x*(1.-lam) - (includeEikonal ?  2.*beta*lam : 0. ));
	s4 =-x*rE*dm/ephig*(1.-lam)*(1.+ctheta)/sqrt(1.+beta);
      }
      else if(fhel==1) {
	s1 =-x*rE*dm*ephig*(1.+lam)*(1.+ctheta)/sqrt(1.+beta);
	s2 =-rE*dm/sqrt(1.+beta)*stheta*
	  ( x*(1.+lam) + (includeEikonal ?  2.*beta*lam : 0. ));
	s3 =-x*rE*dEf*sqrt(1.+beta)*ephig*(1.-lam)*sqr(stheta)/(1.+ctheta);
	s4 = rE*dEf*sqrt(1.+beta)*stheta*
	  ( x*(1.-lam) + (includeEikonal ?  2.*beta*lam : 0. ));
      }
    }
    // anti fermion aling -z
    else {
      if(fhel==0) {
	s1 = -x*rE*dEf*sqrt(1.+beta)/ephig*(1.+lam)*sqr(stheta)/(1.+ctheta);
	s2 = rE*dEf*sqrt(1.+beta)*stheta*
	  ( x*(1.+lam) - (includeEikonal ? 2.*beta*lam : 0. ));
	s3 = x*rE*dm/ephig*(-1.+lam)*(1.+ctheta)/sqrt(1.+beta);
	s4 =-rE*dm/sqrt(1.+beta)*stheta*
	  (+x*(1.-lam) - (includeEikonal ? 2.*beta*lam : 0. ));
      }
      else if(fhel==1) {
	s1 =-rE*dm/sqrt(1.+beta)*stheta*
	  ( x*(1.+lam) + (includeEikonal ? 2.*beta*lam : 0. ));
	s2 = x*rE*dm*ephig*(lam+1.)*(1.+ctheta)/sqrt(1.+beta);
	s3 = rE*dEf*sqrt(1.+beta)*stheta*
	  ( x*(1.-lam) + (includeEikonal ? 2.*beta*lam : 0. ));
	s4 =-x*rE*dEf*sqrt(1.+beta)*ephig*(lam-1.)*sqr(stheta)/(1.+ctheta);
      }
    }
  }
  s1 *= -fact;
  s2 *= -fact;
  s3 *= -fact;
  s4 *= -fact; 
  return SpinorBarWaveFunction(pout,out,s1,s2,s3,s4);
}
