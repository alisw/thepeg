// -*- C++ -*-
//
// VSSVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VSSVertex class.
//

#include "VSSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;
    
// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<VSSVertex,GeneralVSSVertex>
describeThePEGVSSVertex("ThePEG::VSSVertex", "libThePEG.so");

void VSSVertex::Init() {
      
static ClassDocumentation<VSSVertex> documentation
  ("The VSSVertex class is hte implementation of the"
   "vector-scalar-scalar vertex for helicity amplitude calculations."
   " all such vertices should inherit from it"); 
}
