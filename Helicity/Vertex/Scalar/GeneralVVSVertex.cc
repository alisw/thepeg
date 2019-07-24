// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralVVSVertex class.
//

#include "GeneralVVSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<GeneralVVSVertex,AbstractVVSVertex>
describeThePEGGeneralVVSVertex("ThePEG::GeneralVVSVertex", "libThePEG.so");

void GeneralVVSVertex::Init() {

  static ClassDocumentation<GeneralVVSVertex> documentation
    ("The GeneralVVSVertex class implements a general form of"
     " the vector-vector-scalar interaction");

}

Complex GeneralVVSVertex::evaluate(Energy2 q2,const VectorWaveFunction & vec1,
				   const VectorWaveFunction & vec2,
				   const ScalarWaveFunction & sca) {
  Lorentz5Momentum pSca = sca.momentum();
  Lorentz5Momentum pvec1 = vec1.momentum();
  Lorentz5Momentum pvec2 = vec2.momentum();
  // calculate kinematics
  if(kinematics()) calculateKinematics(pSca,pvec1,pvec2);
  // calculate coupling
  setCoupling(q2, vec1.particle(), vec2.particle(), sca.particle());
  Complex e1e2(vec1.wave().dot(vec2.wave()));
  complex<Energy> e1p1(vec1.wave().dot(pvec1));
  complex<Energy> e1p2(vec1.wave().dot(pvec2));
  complex<Energy> e2p1(vec2.wave().dot(pvec1));
  complex<Energy> e2p2(vec2.wave().dot(pvec2));
  complex<Energy2> p1p2(invariant(1,2));
  LorentzPolarizationVectorE eps = epsilon(vec1.wave(),vec2.wave(),pvec2);
  complex<Energy2> p1Ep2 = eps.dot(pvec1);

  Complex output = UnitRemoval::InvE2 * (_a00*e1e2*p1p2 + _aEp*p1Ep2 + 
  _a11*e1p1*e2p1 + _a12*e1p1*e2p2 + _a21*e1p2*e2p1+ _a22*e1p2*e2p2);
  return -norm()*Complex(0.,1.) * sca.wave() * output;
}

ScalarWaveFunction GeneralVVSVertex::evaluate(Energy2 q2,int iopt, tcPDPtr out,
					      const VectorWaveFunction & vec1,
					      const VectorWaveFunction & vec2,
					      complex<Energy> mass,
					      complex<Energy> width) {
  // pointers to the particle data objects
  tcPDPtr Pvec1(vec1.particle());
  tcPDPtr Pvec2(vec2.particle());
  Lorentz5Momentum pvec1 = vec1.momentum();
  Lorentz5Momentum pvec2 = vec2.momentum();
  Lorentz5Momentum pout = pvec1 + pvec2;
  pout.rescaleMass();
  // calculate kinematics if needed
  if(kinematics()) calculateKinematics(-pout,pvec1,pvec2);
  // calculate coupling
  setCoupling(q2,Pvec1,Pvec2,out);
  // propagator
  Complex prop = propagator(iopt,pout.m2(),out,mass,width);
  // lorentz part
  Complex e1e2(vec1.wave().dot(vec2.wave()));
  complex<Energy> e1p1(vec1.wave().dot(pvec1));
  complex<Energy> e1p2(vec1.wave().dot(pvec2));
  complex<Energy> e2p1(vec2.wave().dot(pvec1));
  complex<Energy> e2p2(vec2.wave().dot(pvec2));
  complex<Energy2> p1p2(invariant(1,2));
  LorentzPolarizationVectorE eps = epsilon(vec1.wave(),vec2.wave(),pvec2);
  complex<Energy2> p1Ep2 = eps.dot(pvec1);
  Complex output = UnitRemoval::InvE2 * (_a00*e1e2*p1p2 + _aEp*p1Ep2 + 
  _a11*e1p1*e2p1 + _a12*e1p1*e2p2 + _a21*e1p2*e2p1+ _a22*e1p2*e2p2);
  output *= norm()*prop;
  return ScalarWaveFunction(pout,out,output);
}

VectorWaveFunction GeneralVVSVertex::evaluate(Energy2 q2,int iopt,tcPDPtr out,
					      const VectorWaveFunction & vec,
					      const ScalarWaveFunction & sca,
					      complex<Energy> mass,
					      complex<Energy> width) {
  Lorentz5Momentum pSca = sca.momentum();
  Lorentz5Momentum pvec1 = vec.momentum()+sca.momentum();
  Lorentz5Momentum pvec2 = vec.momentum();
  // calculate kinematics
  if(kinematics()) calculateKinematics(pSca,-pvec1,pvec2);
  // calculate coupling
  setCoupling(q2, out, vec.particle(), sca.particle());
  // prefactor
  Energy2 p2    = pvec1.m2();
  if(mass.real() < ZERO) mass   = out->mass();
  complex<Energy2> mass2 = sqr(mass);
  Complex fact = -norm()* sca.wave() * propagator(iopt,p2,out,mass,width);
  // vertex as polarization vector
  complex<Energy> e2p1(-vec.wave().dot(pvec1));
  complex<Energy> e2p2(vec.wave().dot(pvec2));
  complex<Energy2> p1p2(invariant(1,2));
  LorentzPolarizationVector pv =  (LorentzPolarizationVector(UnitRemoval::InvE2*_a00*p1p2*vec.wave()) -
				   LorentzPolarizationVector(UnitRemoval::InvE2*_a11*e2p1*pvec1) -
				   LorentzPolarizationVector(UnitRemoval::InvE2*_a12*e2p2*pvec1) +
				   LorentzPolarizationVector(UnitRemoval::InvE2*_a21*e2p1*pvec2) +
				   LorentzPolarizationVector(UnitRemoval::InvE2*_a22*e2p2*pvec2) +
				   LorentzPolarizationVector(UnitRemoval::InvE2*_aEp*epsilon(pvec1,vec.wave(),pvec2)));
  // evaluate the wavefunction
  LorentzPolarizationVector vect;
  // massless case
  if(mass.real()==ZERO) {
    vect = fact*pv;
  }
  // massive case
  else {
    complex<InvEnergy> dot = pv.dot(pvec1)/mass2;
    vect = fact*(pv-dot*pvec1);
  }
  return VectorWaveFunction(pvec1,out,vect);
}

