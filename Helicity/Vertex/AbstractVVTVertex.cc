// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVTVertex class.
//

#include "AbstractVVTVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractVVTVertex,VertexBase>
describeThePEGAbstractVVTVertex("ThePEG::AbstractVVTVertex", "libThePEG.so");

void AbstractVVTVertex::Init() {

  static ClassDocumentation<AbstractVVTVertex> documentation
    ("The AbstractVVTVertex class is the base class for all "
     "vector-vector-tensor interactions in ThEPEG");

}

