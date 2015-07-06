// -*- C++ -*-
//
// FermionSpinInfo.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FermionSpinInfo class.
//
// Author: Peter Richardson
//

#include "FermionSpinInfo.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

NoPIOClassDescription<FermionSpinInfo> FermionSpinInfo::initFermionSpinInfo;
// Definition of the static class description member.

void FermionSpinInfo::Init() {}

void FermionSpinInfo::transform(const LorentzMomentum & m,
				const LorentzRotation & r) {
  if(isNear(m)) {
    for(unsigned int ix=0;ix<2;++ix) _currentstates[ix].transform(r);
    SpinInfo::transform(m,r);
  }
}

EIPtr FermionSpinInfo::clone() const {
  tcSpinPtr temp=this;
  return const_ptr_cast<SpinPtr>(temp);
}
