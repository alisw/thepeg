// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralVSSVertex class.
//

#include "GeneralVSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<GeneralVSSVertex,AbstractVSSVertex>
describeHelicityGeneralVSSVertex("ThePEG::GeneralVSSVertex", "libThePEG.so");

void GeneralVSSVertex::Init() {

  static ClassDocumentation<GeneralVSSVertex> documentation
    ("The GeneralVSSVertex class implements a general form of the vector-scalar-scalar interaction");

}

// evaluate the vertex
Complex GeneralVSSVertex::evaluate(Energy2 q2, const VectorWaveFunction & vec,
				   const ScalarWaveFunction & sca1,
				   const ScalarWaveFunction & sca2) {
  // calculate the coupling
  setCoupling(q2,vec.particle(),sca1.particle(),sca2.particle());
  // calculate the vertex
  return UnitRemoval::InvE * -Complex(0.,1.) * norm() * sca1.wave()*sca2.wave()*
    vec.wave().dot(a_*sca1.momentum()+b_*sca2.momentum());
}

// off-shell vector
VectorWaveFunction GeneralVSSVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
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
  LorentzPolarizationVector vec =
    -UnitRemoval::InvE * fact * (b_*sca2.momentum()+a_*sca1.momentum());
  // massive outgoing vector
  if(mass.real()!=ZERO) {
    // first the dot product for the second term
    complex<InvEnergy> dot = vec*pout/mass2;
    vec -= dot*pout;
  }
  return VectorWaveFunction(pout,out,vec);
}

// return an off-shell scalar
ScalarWaveFunction GeneralVSSVertex::evaluate(Energy2 q2, int iopt, tcPDPtr out,
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
  fact = Complex(UnitRemoval::InvE * fact*vec.wave().dot(a_*sca.momentum()-b_*pout));
  return ScalarWaveFunction(pout,out,fact);
}
