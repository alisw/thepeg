// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractSSSTVertex class.
//

#include "AbstractSSSTVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractSSSTVertex,VertexBase>
describeThePEGAbstractSSSTVertex("ThePEG::AbstractSSSTVertex", "libThePEG.so");

void AbstractSSSTVertex::Init() {

  static ClassDocumentation<AbstractSSSTVertex> documentation
    ("The AbstractSSSTVertex class is the base class for scalar-scalar-scalar-tensor"
     " interactions in ThePEG");

}

