// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralVVSVertex class.
//

#include "GeneralVVSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<GeneralVVSVertex> GeneralVVSVertex::initGeneralVVSVertex;
// Definition of the static class description member.

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
  if(kinematics()) calculateKinematics(pout,pvec1,pvec2);
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
  output *=  -norm()*prop;
  return ScalarWaveFunction(pout,out,output);
}

VectorWaveFunction GeneralVVSVertex::evaluate(Energy2 ,int ,tcPDPtr ,
					      const VectorWaveFunction & ,
					      const ScalarWaveFunction & ,
					      complex<Energy>,
					      complex<Energy>) {
  throw Exception() << "GeneralVVSVertex::evaluate() only implemented for the "
		    << "member which returns the amplitude, "
		    << "and off-shell scalar not the off-shell vector"
		    << Exception::runerror;
}

