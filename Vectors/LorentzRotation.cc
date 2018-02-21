// -*- C++ -*-
//
// LorentzRotation.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LorentzRotation class.
//

#include "LorentzRotation.h"

using namespace ThePEG;

// output operator
std::ostream & LorentzRotation::print( std::ostream & os ) const
{
  os <<    "Spin 1   Transform: \n " << _one 
     << "\n Spin 1/2 Transform: \n " << _half << "\n";
  return os;
}
