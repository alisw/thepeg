// -*- C++ -*-
//
// RFVVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RFVVertex class.
//

#include "RFVVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// Definition of the static class description member
// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<RFVVertex,AbstractRFVVertex>
describeThePEGRFVVertex("ThePEG::RFVVertex", "libThePEG.so");
    
void RFVVertex::Init() {
      
  static ClassDocumentation<RFVVertex> documentation
    ("The RFVVertex class implements the helicity amplitude"
     "calculations for a spin-3/2 fermion-fantifermion gauge boson vertex. Any   "
     "implementation of such a vertex should inherit from in and implement"
     " the virtual setCoupling member to calculate the coupling");
}


Complex RFVVertex::evaluate(Energy2 q2,const RSSpinorWaveFunction & sp,
			    const SpinorBarWaveFunction & sbar,
			    const VectorWaveFunction & vec) {
  // calculate the couplings
  setCoupling(q2,sp.particle(),sbar.particle(),vec.particle());
  LorentzSpinor<double> wdot1 = sp.wave().dot(vec.wave());
  Complex lS1 = wdot1. leftScalar(sbar.wave());
  Complex rS1 = wdot1.rightScalar(sbar.wave());
  LorentzSpinor<double> wdot2 = sp.wave().dot(sbar.momentum());
  Complex lS2 = wdot2. leftCurrent(sbar.wave()).dot(vec.wave());
  Complex rS2 = wdot2.rightCurrent(sbar.wave()).dot(vec.wave());
  swap(lS2,rS2);
  Complex dot = sbar.momentum().dot(vec.wave())*UnitRemoval::InvE; 
  Complex lS3 = wdot2. leftScalar(sbar.wave())*dot;
  Complex rS3 = wdot2.rightScalar(sbar.wave())*dot;
  return Complex(0.,1.)*norm()*
    (lS1*left()[0]+rS1*right()[0]+
     lS2*left()[1]+rS2*right()[1]+
     lS3*left()[2]+rS3*right()[2]);
}

Complex RFVVertex::evaluate(Energy2 q2,const SpinorWaveFunction & sp,
			    const RSSpinorBarWaveFunction & sbar,
			    const VectorWaveFunction & vec) {
  // calculate the couplings
  setCoupling(q2,sbar.particle(),sp.particle(),vec.particle());
  LorentzSpinorBar<double> wdot1 = sbar.wave().dot(vec.wave());
  Complex lS1 = sp.wave(). leftScalar(wdot1);
  Complex rS1 = sp.wave().rightScalar(wdot1);
  LorentzSpinorBar<double> wdot2 = sbar.wave().dot(sp.momentum());
  Complex lS2 = sp.wave(). leftCurrent(wdot2).dot(vec.wave());
  Complex rS2 = sp.wave().rightCurrent(wdot2).dot(vec.wave());
  Complex dot = sbar.momentum().dot(vec.wave())*UnitRemoval::InvE; 
  Complex lS3 = sp.wave(). leftScalar(wdot2)*dot;
  Complex rS3 = sp.wave().rightScalar(wdot2)*dot;
  return Complex(0.,1.)*norm()*
    (lS1*left()[0]+rS1*right()[0]+
     lS2*left()[1]+rS2*right()[1]+
     lS3*left()[2]+rS3*right()[2]);
}

SpinorBarWaveFunction RFVVertex::evaluate(Energy2 ,int ,tcPDPtr ,
					  const RSSpinorBarWaveFunction & ,
					  const VectorWaveFunction & ,
					  complex<Energy> , complex<Energy> ) {
  assert(false);
  return SpinorBarWaveFunction();
}

RSSpinorBarWaveFunction RFVVertex::evaluate(Energy2 ,int ,tcPDPtr ,
					    const SpinorBarWaveFunction & ,
					    const VectorWaveFunction & ,
					    complex<Energy> , complex<Energy> ) {
  assert(false);
  return RSSpinorBarWaveFunction();
}

VectorWaveFunction RFVVertex::evaluate(Energy2 ,int ,tcPDPtr ,
				       const RSSpinorWaveFunction & ,
				       const SpinorBarWaveFunction & ,
				       complex<Energy> , complex<Energy> ) {
  assert(false);
  return VectorWaveFunction();
}

VectorWaveFunction RFVVertex::evaluate(Energy2 ,int ,tcPDPtr ,
				       const SpinorWaveFunction & ,
				       const RSSpinorBarWaveFunction & ,
				       complex<Energy> , complex<Energy> ) {
  assert(false);
  return VectorWaveFunction();
}

RSSpinorWaveFunction RFVVertex::evaluate(Energy2 ,int ,tcPDPtr ,
					 const SpinorWaveFunction & ,
					 const VectorWaveFunction & ,
					 complex<Energy> , complex<Energy> ) {
  assert(false);
  return RSSpinorWaveFunction();
}

SpinorWaveFunction RFVVertex::evaluate(Energy2 ,int ,tcPDPtr ,
				       const RSSpinorWaveFunction & ,
				       const VectorWaveFunction & ,
				       complex<Energy> , complex<Energy> ) {
  assert(false);
  return SpinorWaveFunction();
}
