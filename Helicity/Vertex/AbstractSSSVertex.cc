// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractSSSVertex class.
//

#include "AbstractSSSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractSSSVertex,VertexBase>
describeThePEGAbstractSSSVertex("ThePEG::AbstractSSSVertex", "libThePEG.so");

void AbstractSSSVertex::Init() {

  static ClassDocumentation<AbstractSSSVertex> documentation
    ("The AbstractSSSVertex class is the base class for all "
     "scalar-scalar-scalar interactions");

}

