// -*- C++ -*-
//
// GeneralFFVVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralFFVVertex class.
//

#include "GeneralFFVVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// Definition of the static class description member
// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<GeneralFFVVertex,AbstractFFVVertex>
describeThePEGGeneralFFVVertex("ThePEG::GeneralFFVVertex", "libThePEG.so");
    
void GeneralFFVVertex::Init() {
      
  static ClassDocumentation<GeneralFFVVertex> documentation
    ("The GeneralFFVVertex class implements the helicity amplitude"
     "calculations for a fermion-fantifermion gauge boson vertex. Any   "
     "implementation of such a vertex should inherit from in and implement"
     " the virtual setCoupling member to calculate the coupling");

}

// evalulate the full vertex
Complex GeneralFFVVertex::evaluate(Energy2 q2,
				   const SpinorWaveFunction & sp, 
				   const SpinorBarWaveFunction & sbar,
				   const VectorWaveFunction & vec) {
  // first calculate the couplings
  setCoupling(q2,sp.particle(),sbar.particle(),vec.particle());
  const Complex ii(0.,1.);
  const complex<InvEnergy> zero(ZERO,ZERO);
  // useful combinations of the polarization vector components
  Complex e0p3=vec.t()+vec.z();
  Complex e0m3=vec.t()-vec.z();
  Complex e1p2=vec.x()+ii*vec.y();
  Complex e1m2=vec.x()-ii*vec.y();
  complex<Energy> p0p3=vec. e()+vec.pz();
  complex<Energy> p0m3=vec. e()-vec.pz();
  complex<Energy> p1p2=vec.px()+ii*vec.py();
  complex<Energy> p1m2=vec.px()-ii*vec.py();
  // get both the spinors in the same representation (using the default choice)
  Complex vertex(0.);
  // first the left piece as this is virtually always needed
  if(_left!=0.) {
    vertex = _left*(+sbar.s3()*(sp.s1()*e0p3+sp.s2()*e1m2)
		    +sbar.s4()*(sp.s1()*e1p2+sp.s2()*e0m3));
  }
  // then the right piece (often not needed eg W vertex)
  if(_right!=0.) {
    vertex += _right*(+sbar.s1()*(sp.s3()*e0m3-sp.s4()*e1m2)
		      -sbar.s2()*(sp.s3()*e1p2-sp.s4()*e0p3));
  }
  // left sigma piece
  if(_leftSigma!=zero) {
    vertex -= Complex(ii * _leftSigma *
		       (sbar.s1()*sp.s1()*(-vec. e()*vec.z() + vec.pz()*vec.t() 
					   -ii*(vec.py()*vec.x()-vec.px()*vec.y()))+
			sbar.s2()*sp.s1()*(   p1p2*e0p3 - p0p3*e1p2 ) +
			sbar.s1()*sp.s2()*( + p1m2*e0m3 - p0m3*e1m2 ) +
			sbar.s2()*sp.s2()*( vec. e()*vec.z() - vec.pz()*vec.t()
					    -ii*(vec.px()*vec.y()-vec.py()*vec.x()))));
  }
  // right sigma piece
  if(_rightSigma!=zero) {
    vertex += Complex(ii * _rightSigma *
		      (sbar.s3()*sp.s3()*(-vec.e()*vec.z()+vec.pz()*vec.t()
					  -ii*(vec.px()*vec.y()-vec.py()*vec.x()))+
		       sbar.s4()*sp.s3()*(   p1p2*e0m3 - p0m3*e1p2 ) +
		       sbar.s3()*sp.s4()*(   p1m2*e0p3 - p0p3*e1m2 ) +
		       sbar.s4()*sp.s4()*( vec.e()*vec.z()-vec.pz()*vec.t()
					   +ii*(vec.px()*vec.y()-vec.py()*vec.x()))));
  }
  vertex *= ii;
  // final factors
  return vertex*norm();
}

// evaluate an off-shell spinor
SpinorWaveFunction GeneralFFVVertex::evaluate(Energy2 q2, int iopt,tcPDPtr  out,
					      const SpinorWaveFunction & sp,
					      const VectorWaveFunction &vec,
					      complex<Energy> mass,
					      complex<Energy> width) {
  // extract the pointers to the particle data objects
  tcPDPtr  Psp=sp.particle();
  tcPDPtr  Pvec=vec.particle();
  // first calculate the couplings
  setCoupling(q2,Psp,out,Pvec);
  const Complex ii(0.,1.);
  const complex<InvEnergy> zero(ZERO,ZERO);
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = sp.momentum()+vec.momentum();
  // now evaluate the contribution
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
  complex<Energy> p0p3 = pout.e() +    pout.z();
  complex<Energy> p0m3 = pout.e() -    pout.z();
  complex<Energy> p1p2 = pout.x() + ii*pout.y();
  complex<Energy> p1m2 = pout.x() - ii*pout.y();
  // complex nos for for the spinor
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  Complex a1 = fact*( sp.s3()*e0m3-sp.s4()*e1m2);
  Complex a2 = fact*(-sp.s3()*e1p2+sp.s4()*e0p3);
  Complex a3 = fact*( sp.s1()*e0p3+sp.s2()*e1m2);
  Complex a4 = fact*( sp.s1()*e1p2+sp.s2()*e0m3);
  // left piece
  if(_left!=0.) {
    s1 += Complex(UnitRemoval::InvE * _left * (p0m3*a3-p1m2*a4));
    s2 += Complex(UnitRemoval::InvE * _left * (-p1p2*a3+p0p3*a4));
    s3 += Complex(UnitRemoval::InvE * _left * a3*mass);
    s4 += Complex(UnitRemoval::InvE * _left * a4*mass);
  }
  // right piece
  if(_right!=0.) {
    s1 += Complex(UnitRemoval::InvE * _right * a1*mass);
    s2 += Complex(UnitRemoval::InvE * _right * a2*mass);
    s3 += Complex(UnitRemoval::InvE * _right * (p0p3*a1+p1m2*a2));
    s4 += Complex(UnitRemoval::InvE * _right * (p1p2*a1+p0m3*a2));
  }
  complex<Energy> p0p3b = vec. e() +    vec.pz();
  complex<Energy> p0m3b = vec. e() -    vec.pz();
  complex<Energy> p1p2b = vec.px() + ii*vec.py();
  complex<Energy> p1m2b = vec.px() - ii*vec.py();
  complex<Energy> b1 = fact*( sp.s3()*p0m3b - sp.s4()*p1m2b);
  complex<Energy> b2 = fact*(-sp.s3()*p1p2b + sp.s4()*p0p3b);
  complex<Energy> b3 = fact*( sp.s1()*p0p3b + sp.s2()*p1m2b);
  complex<Energy> b4 = fact*( sp.s1()*p1p2b + sp.s2()*p0m3b);
  // left sigma piece
  if(_leftSigma!=zero) {
    s1 += Complex(-0.5*ii*UnitRemoval::InvE *mass*_leftSigma*
		  ( - a3*p0m3b + a4*p1m2b + b3*e0m3 - b4*e1m2));
    s2 += Complex( 0.5*ii*UnitRemoval::InvE *mass*_leftSigma*
		   ( - a3*p1p2b + a4*p0p3b + b3*e1p2 - b4*e0p3));
    s3 += Complex(-0.5*ii*UnitRemoval::InvE*_leftSigma*
		  ( + p0p3*( - a3*p0m3b + a4*p1m2b + b3*e0m3 - b4*e1m2 )
		    + p1m2*(   a3*p1p2b - a4*p0p3b - b3*e1p2 + b4*e0p3 ) ));
    s4 += Complex( 0.5*ii*UnitRemoval::InvE*_leftSigma*
		   ( + p1p2*( + a3*p0m3b - a4*p1m2b - b3*e0m3 + b4*e1m2 )
		     + p0m3*( - a3*p1p2b + a4*p0p3b + b3*e1p2 - b4*e0p3 ) ));
  }
  // right sigma piece
  if(_rightSigma!=zero) {
    s1 +=Complex(  0.5*ii*UnitRemoval::InvE *_rightSigma*
		   ( + p0m3*( + a1*p0p3b + a2*p1m2b - b1*e0p3 - b2*e1m2 )
		     + p1m2*( - a1*p1p2b - a2*p0m3b + b1*e1p2 + b2*e0m3 ) ));
    s2 +=Complex( -0.5*ii*UnitRemoval::InvE *_rightSigma*
		  ( + p1p2*( + a1*p0p3b + a2*p1m2b - b1*e0p3 - b2*e1m2 )
		    + p0p3*( - a1*p1p2b - a2*p0m3b + b1*e1p2 + b2*e0m3 ) ));
    s3 += Complex( 0.5*ii*UnitRemoval::InvE *mass*_rightSigma*
		   (  a1*p0p3b + a2*p1m2b - b1*e0p3 - b2*e1m2 ));
    s4 +=Complex( -0.5*ii*UnitRemoval::InvE *mass*_rightSigma*
		  ( -a1*p1p2b - a2*p0m3b + b2*e0m3 + b1*e1p2 ));
  }
  // return the wavefunction
  return SpinorWaveFunction(pout,out,s1,s2,s3,s4);
}
  
// evaluate an off-shell SpinorBar
SpinorBarWaveFunction GeneralFFVVertex::evaluate(Energy2 q2,int iopt,tcPDPtr  out,
						 const SpinorBarWaveFunction & sbar,
						 const VectorWaveFunction & vec,
						 complex<Energy> mass,
						 complex<Energy> width) {
  // first calculate the couplings
  setCoupling(q2,out,sbar.particle(),vec.particle());
  const Complex ii(0.,1.);
  const complex<InvEnergy> zero(ZERO,ZERO);
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = sbar.momentum()+vec.momentum();
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
  complex<Energy> p1p2=pout.x() + ii*pout.y();
  complex<Energy> p1m2=pout.x() - ii*pout.y();
  complex<Energy> p0p3=pout.e() +    pout.z();
  complex<Energy> p0m3=pout.e() -    pout.z();
  // complex numbers for the spinor
  Complex s1(0.),s2(0.),s3(0.),s4(0.);
  Complex a1 = fact*( sbar.s3()*e0p3+sbar.s4()*e1p2);
  Complex a2 = fact*( sbar.s3()*e1m2+sbar.s4()*e0m3);
  Complex a3 = fact*( sbar.s1()*e0m3-sbar.s2()*e1p2);
  Complex a4 = fact*(-sbar.s1()*e1m2+sbar.s2()*e0p3);
  complex<Energy> p0p3b = vec. e() +    vec.pz();
  complex<Energy> p0m3b = vec. e() -    vec.pz();
  complex<Energy> p1p2b = vec.px() + ii*vec.py();
  complex<Energy> p1m2b = vec.px() - ii*vec.py();
  complex<Energy> b1 = fact*( sbar.s3()*p0p3b+sbar.s4()*p1p2b);
  complex<Energy> b2 = fact*( sbar.s3()*p1m2b+sbar.s4()*p0m3b);
  complex<Energy> b3 = fact*( sbar.s1()*p0m3b-sbar.s2()*p1p2b);
  complex<Energy> b4 = fact*(-sbar.s1()*p1m2b+sbar.s2()*p0p3b);
  // left piece
  if(_left!=0.) {
    s1 += Complex(UnitRemoval::InvE*_left*a1*mass);
    s2 += Complex(UnitRemoval::InvE*_left*a2*mass);
    s3 += Complex(UnitRemoval::InvE*_left*(-p0m3*a1+p1p2*a2));
    s4 += Complex(UnitRemoval::InvE*_left*(+p1m2*a1-p0p3*a2));
  }
  // right piece
  if(_right!=0.) {
    s1 += Complex(UnitRemoval::InvE*_right*(-p0p3*a3-p1p2*a4));
    s2 += Complex(UnitRemoval::InvE*_right*(-p1m2*a3-p0m3*a4));
    s3 += Complex(UnitRemoval::InvE*_right*a3*mass);
    s4 += Complex(UnitRemoval::InvE*_right*a4*mass); 
  }
  // left sigma piece
  if(_leftSigma!=zero) {
    s1 +=Complex(  0.5*ii*UnitRemoval::InvE*_leftSigma*mass*
		   ( + a3*p0p3b + a4*p1p2b - b3*e0p3 - b4*e1p2));
    s2 +=Complex( -0.5*ii*UnitRemoval::InvE*_leftSigma*mass*
		  ( - a3*p1m2b - a4*p0m3b + b3*e1m2 + b4*e0m3));
    s3 +=Complex( -0.5*ii*UnitRemoval::InvE*_leftSigma*
		  (+p0m3*( + a3*p0p3b + a4*p1p2b - b3*e0p3 - b4*e1p2)
		   +p1p2*( - a3*p1m2b - a4*p0m3b + b3*e1m2 + b4*e0m3)));
    s4 +=Complex( +0.5*ii*UnitRemoval::InvE*_leftSigma*
		  (+p1m2*( + a3*p0p3b + a4*p1p2b - b3*e0p3 - b4*e1p2)
		   +p0p3*( - a3*p1m2b - a4*p0m3b + b3*e1m2 + b4*e0m3)));
  }
  // right sigma piece
  if(_rightSigma!=zero) {
    s1 +=Complex( +0.5*ii*UnitRemoval::InvE*_rightSigma*
		  (+p0p3*( - a1*p0m3b + a2*p1p2b + b1*e0m3 - b2*e1p2)
		   +p1p2*( + a1*p1m2b - a2*p0p3b - b1*e1m2 + b2*e0p3)));
    s2 +=Complex( -0.5*ii*UnitRemoval::InvE*_rightSigma*
		  (+p1m2*( + a1*p0m3b - a2*p1p2b - b1*e0m3 + b2*e1p2)
		   +p0m3*( - a1*p1m2b + a2*p0p3b + b1*e1m2 - b2*e0p3)));
    s3 +=Complex( -0.5*ii*UnitRemoval::InvE*_rightSigma*mass*
		  ( - a1*p0m3b + a2*p1p2b + b1*e0m3 - b2*e1p2 ));
    s4 +=Complex(  0.5*ii*UnitRemoval::InvE*_rightSigma*mass*
		   ( - a1*p1m2b + a2*p0p3b + b1*e1m2 - b2*e0p3 ));
  }
  return SpinorBarWaveFunction(pout,out,s1,s2,s3,s4);
}

// off-shell vector
VectorWaveFunction GeneralFFVVertex::evaluate(Energy2 q2,int iopt,tcPDPtr  out,
					      const SpinorWaveFunction & sp,
					      const SpinorBarWaveFunction & sbar,
					      complex<Energy> mass,
					      complex<Energy> width) {
  // first calculate the couplings
  setCoupling(q2,sp.particle(),sbar.particle(),out);
  const Complex ii(0.,1.);
  const complex<InvEnergy> zero(ZERO,ZERO);
  // work out the momentum of the off-shell particle
  Lorentz5Momentum pout = sbar.momentum()+sp.momentum();
  // overall factor
  Energy2 p2 = pout.m2();
  Complex fact = normPropagator(iopt,p2,out,mass,width);
  // momentum components
  if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
  complex<Energy2> mass2 = sqr(mass);
  // the vector for the fermion-antifermion
  Complex vec[4];
  complex<Energy> p1p2=pout.x() + ii*pout.y();
  complex<Energy> p1m2=pout.x() - ii*pout.y();
  complex<Energy> p0p3=pout.e() +    pout.z();
  complex<Energy> p0m3=pout.e() -    pout.z();
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
  // left sigma
  if(_leftSigma==zero) {
    vec[0] +=Complex( -0.5*ii*_leftSigma*
		      (+sp.s1()*sbar.s1()*(p1p2-p1m2)+2.*sp.s1()*sbar.s2()*p0p3
		       -sp.s2()*sbar.s2()*(p1p2-p1m2)+2.*sp.s2()*sbar.s1()*p0m3));
    vec[1] +=Complex( -0.5   *_leftSigma*
		      (+sp.s1()*sbar.s1()*(p1m2-p1p2)+2.*sp.s1()*sbar.s2()*p0p3
		       +sp.s2()*sbar.s2()*(p1p2+p1m2)-2.*sp.s2()*sbar.s1()*p0m3));
    vec[2] +=Complex( -0.5*ii*_leftSigma*
		      (+sp.s1()*sbar.s1()*(p0p3+p0m3)-2.*sp.s1()*sbar.s2()*p1p2
		       -sp.s2()*sbar.s2()*(p0p3+p0m3)+2.*sp.s2()*sbar.s1()*p1m2));
    vec[3] +=Complex(  0.5*ii*_leftSigma*
		       (-sp.s1()*sbar.s1()*(p0p3-p0m3)-2.*sp.s1()*sbar.s2()*p1p2
			+sp.s2()*sbar.s2()*(p0p3-p0m3)-2.*sp.s2()*sbar.s1()*p1m2));
  }
  // right sigma
  if(_rightSigma==zero) {
    vec[0] +=Complex(  0.5*ii*_rightSigma*
		       (-sp.s3()*sbar.s3()*(p1p2-p1m2)+2.*sp.s3()*sbar.s4()*p0m3
			+sp.s4()*sbar.s4()*(p1p2-p1m2)+2.*sp.s4()*sbar.s3()*p0p3));
    vec[1] +=Complex(  0.5*   _rightSigma*
		       (-sp.s3()*sbar.s3()*(p1p2+p1m2)-2.*sp.s3()*sbar.s4()*p0m3
			+sp.s4()*sbar.s4()*(p1p2+p1m2)+2.*sp.s4()*sbar.s3()*p0p3));
    vec[2] +=Complex(  0.5*ii*_rightSigma*
		       (+sp.s3()*sbar.s3()*(p0p3+p0m3)+2.*sp.s3()*sbar.s4()*p1p2
			-sp.s4()*sbar.s4()*(p0p3+p0m3)-2.*sp.s4()*sbar.s3()*p1m2));
    vec[3] +=Complex(  -0.5*ii*_rightSigma*
		       (-sp.s3()*sbar.s3()*(p0p3-p0m3)-2.*sp.s3()*sbar.s4()*p1p2
			+sp.s4()*sbar.s4()*(p0p3-p0m3)-2.*sp.s4()*sbar.s3()*p1m2));
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
