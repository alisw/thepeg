// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractSSTVertex class.
//

#include "AbstractSSTVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractSSTVertex,VertexBase>
describeThePEGAbstractSSTVertex("ThePEG::AbstractSSTVertex", "libThePEG.so");

void AbstractSSTVertex::Init() {

  static ClassDocumentation<AbstractSSTVertex> documentation
    ("The AbstractSSTVertex class is the base class for scalar-scalar-tensor"
     "interactions in ThePEG.");

}

