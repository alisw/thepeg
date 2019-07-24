// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVSSVertex class.
//

#include "AbstractVVSSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractVVSSVertex,VertexBase>
describeThePEGAbstractVVSSVertex("ThePEG::AbstractVVSSVertex", "libThePEG.so");

void AbstractVVSSVertex::Init() {

  static ClassDocumentation<AbstractVVSSVertex> documentation
    ("The AbstractVVSSVertex class is the base class for "
     "vector-vector-scalar-scalar interactions in ThePEG");

}

