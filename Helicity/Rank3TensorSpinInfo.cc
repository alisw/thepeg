// -*- C++ -*-
//
// Rank3TensorSpinInfo.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Rank3TensorSpinInfo class.
//
// Author: Peter Richardson
//

#include "Rank3TensorSpinInfo.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<Rank3TensorSpinInfo,SpinInfo>
describeThePEGRank3TensorSpinInfo("ThePEG::Rank3TensorSpinInfo", "libThePEG.so");

void Rank3TensorSpinInfo::Init() {}

void Rank3TensorSpinInfo::transform(const LorentzMomentum & m,
				    const LorentzRotation & r) {
  if(isNear(m)) {
    for(unsigned int ix=0;ix<7;++ix) _currentstates[ix].transform(r.one());
    SpinInfo::transform(m,r);
  }
}

EIPtr Rank3TensorSpinInfo::clone() const {
  tcSpinPtr temp=this;
  return const_ptr_cast<SpinPtr>(temp);
}
