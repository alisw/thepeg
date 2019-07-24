// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVVTVertex class.
//

#include "AbstractVVVTVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractVVVTVertex,VertexBase>
describeThePEGAbstractVVVTVertex("ThePEG::AbstractVVVTVertex", "libThePEG.so");

void AbstractVVVTVertex::Init() {

  static ClassDocumentation<AbstractVVVTVertex> documentation
    ("The AbstractVVVTVertex class is the base class for vector-vector-vector-tensor"
     " interactions in ThePEG");

}

