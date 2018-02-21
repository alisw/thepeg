// -*- C++ -*-
//
// FFVTVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFVTVertex class.
//

#include "FFVTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<FFVTVertex> FFVTVertex::initFFVTVertex;
// Definition of the static class description member.

void FFVTVertex::Init() {
  
  static ClassDocumentation<FFVTVertex> documentation
    ("The FFVTVertex class is the implementation of the"
     "helicity amplitude calculation of the fermion-antifermion-vector-tensor"
     "vertex. All such vertices should inherit from it.");

}

// function to evaluate the vertex
Complex FFVTVertex::evaluate(Energy2 q2, const SpinorWaveFunction & sp,
			     const SpinorBarWaveFunction & sbar,
			     const VectorWaveFunction & vec,
			     const TensorWaveFunction & ten) {
  // set the couplings
  setCoupling(q2,sp.particle(),sbar.particle(),
	      vec.particle(),ten.particle());
  Complex ii(0.,1.);
  // vector current
  LorentzPolarizationVector as = 
    left_ *sp.wave(). leftCurrent(sbar.wave()) +
    right_*sp.wave().rightCurrent(sbar.wave());
  // trace of the tensor
  Complex trace = ten.wave().trace();
  // dot product
  Complex dotav = as.dot(vec.wave());
  // product with tensor and current
  LorentzPolarizationVector  preDot = ten.wave(). preDot(as);
  LorentzPolarizationVector postDot = ten.wave().postDot(as);
  Complex tenav = preDot.dot(vec.wave())+postDot.dot(vec.wave());
  // return the vertex
  return ii*0.25*norm()*(tenav-2.*trace*dotav);
}

TensorWaveFunction FFVTVertex::evaluate(Energy2,int , tcPDPtr ,
					const SpinorWaveFunction & ,
					const SpinorBarWaveFunction & ,
					const VectorWaveFunction & ,
					complex<Energy>, complex<Energy>) {
  throw Exception() << "FFVTVertex::evaluate() only implemented for the "
		    << "member which returns the amplitude, "
		    << "not the off-shell wavefunctions"
		    << Exception::runerror;
}

VectorWaveFunction FFVTVertex::evaluate(Energy2 ,int , tcPDPtr ,
					const SpinorWaveFunction & ,
					const SpinorBarWaveFunction & , 
					const TensorWaveFunction &  ,
					complex<Energy>, complex<Energy>) {
  throw Exception() << "FFVTVertex::evaluate() only implemented for the "
		    << "member which returns the amplitude, "
		    << "not the off-shell wavefunctions"
		    << Exception::runerror;
}

SpinorWaveFunction FFVTVertex::evaluate(Energy2 ,int , tcPDPtr ,
					const SpinorWaveFunction & ,
					const VectorWaveFunction & ,
					const TensorWaveFunction &  ,
					complex<Energy>, complex<Energy>) {
  throw Exception() << "FFVTVertex::evaluate() only implemented for the "
		    << "member which returns the amplitude, "
		    << "not the off-shell wavefunctions"
		    << Exception::runerror;
}

SpinorBarWaveFunction FFVTVertex::evaluate(Energy2 ,int , tcPDPtr ,
					   const SpinorBarWaveFunction & ,
					   const VectorWaveFunction & ,
					   const TensorWaveFunction &  ,
					   complex<Energy>, complex<Energy>) {
  throw Exception() << "FFVTVertex::evaluate() only implemented for the "
		    << "member which returns the amplitude, "
		    << "not the off-shell wavefunctions"
		    << Exception::runerror;
}
