// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVSTVertex class.
//

#include "AbstractVVSTVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractVVSTVertex,VertexBase>
describeThePEGAbstractVVSTVertex("ThePEG::AbstractVVSTVertex", "libThePEG.so");

void AbstractVVSTVertex::Init() {

  static ClassDocumentation<AbstractVVSTVertex> documentation
    ("The AbstractVVSTVertex class is the base class for vector-vector-scalar-tensor"
     " interactions in ThePEG");

}

