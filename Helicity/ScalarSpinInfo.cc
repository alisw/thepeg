// -*- C++ -*-
//
// ScalarSpinInfo.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarSpinInfo class.
//
// Author: Peter Richardson
//

#include "ScalarSpinInfo.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<ScalarSpinInfo,SpinInfo>
describeThePEGScalarSpinInfo("ThePEG::ScalarSpinInfo", "libThePEG.so");

void ScalarSpinInfo::Init() {}

void ScalarSpinInfo::transform(const LorentzMomentum & m, 
			       const LorentzRotation & r) {
  if(isNear(m))
    SpinInfo::transform(m,r);
}
