// -*- C++ -*-
//
// RFSVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RFSVertex class.
//

#include "RFSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<RFSVertex,AbstractRFSVertex>
describeThePEGRFSVertex("ThePEG::RFSVertex", "libThePEG.so");
    
void RFSVertex::Init() {

  static ClassDocumentation<RFSVertex> documentation
    ("The RFSVertex class is the implementation of the RFS"
     "vertex. All such vertices shoud inherit from it");
  
}

Complex RFSVertex::evaluate(Energy2 q2,const RSSpinorWaveFunction & sp,
			    const SpinorBarWaveFunction & sbar,
			    const ScalarWaveFunction & sca) {
  // calculate the couplings
  setCoupling(q2,sp.particle(),sbar.particle(),sca.particle());
  LorentzSpinor<double> wdot = sp.wave().dot(sbar.momentum());
  Complex lS = wdot. leftScalar(sbar.wave());
  Complex rS = wdot.rightScalar(sbar.wave());
  swap(lS,rS);
  return Complex(0.,1.)*norm()*sca.wave()*(lS*left()+rS*right());
}

Complex RFSVertex::evaluate(Energy2 q2,const SpinorWaveFunction & sp,
		 const RSSpinorBarWaveFunction & sbar,
		 const ScalarWaveFunction & sca) {
  // calculate the couplings
  setCoupling(q2,sbar.particle(),sp.particle(),sca.particle());
  LorentzSpinorBar<double> wdot = sbar.wave().dot(sp.momentum());
  Complex lS = sp.wave(). leftScalar(wdot);
  Complex rS = sp.wave().rightScalar(wdot);
  return Complex(0.,1.)*norm()*sca.wave()*(lS*left()+rS*right());
}


SpinorWaveFunction RFSVertex::evaluate(Energy2 ,int ,tcPDPtr ,
				       const RSSpinorWaveFunction & , 
				       const ScalarWaveFunction & ,
				       complex<Energy> , complex<Energy> ) {
  assert(false);
  return SpinorWaveFunction();
}

RSSpinorWaveFunction RFSVertex::evaluate(Energy2 ,int ,tcPDPtr ,
					 const SpinorWaveFunction & , 
					 const ScalarWaveFunction & ,
					 complex<Energy> , complex<Energy> ) {
  assert(false);
  return RSSpinorWaveFunction();
}

SpinorBarWaveFunction RFSVertex::evaluate(Energy2 ,int ,tcPDPtr ,
					  const RSSpinorBarWaveFunction & ,
					  const ScalarWaveFunction & ,
					  complex<Energy> , complex<Energy> ) {
  assert(false);
  return SpinorBarWaveFunction();
}

RSSpinorBarWaveFunction RFSVertex::evaluate(Energy2 ,int ,tcPDPtr ,
					    const SpinorBarWaveFunction & ,
					    const ScalarWaveFunction & ,
					    complex<Energy> , complex<Energy> ) {
  assert(false);
  return RSSpinorBarWaveFunction();
}

ScalarWaveFunction RFSVertex::evaluate(Energy2 ,int ,tcPDPtr ,
				       const RSSpinorWaveFunction & , 
				       const SpinorBarWaveFunction & ,
				       complex<Energy> , complex<Energy> ) {
  assert(false);
  return ScalarWaveFunction();
}

ScalarWaveFunction RFSVertex::evaluate(Energy2 ,int ,tcPDPtr ,
				       const SpinorWaveFunction & , 
				       const RSSpinorBarWaveFunction & ,
				       complex<Energy> , complex<Energy> ) {
  assert(false);
  return ScalarWaveFunction();
}
