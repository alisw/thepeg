// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVSSVertex class.
//

#include "AbstractVSSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractVSSVertex,VertexBase>
describeThePEGAbstractVSSVertex("ThePEG::AbstractVSSVertex", "libThePEG.so");

void AbstractVSSVertex::Init() {

  static ClassDocumentation<AbstractVSSVertex> documentation
    ("The AbstractVSSVertex class is the base class for all "
     "vector-scalar-scalar interactions in ThePEG.");

}

