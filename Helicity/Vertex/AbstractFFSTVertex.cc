// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractFFSTVertex class.
//

#include "AbstractFFSTVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractFFSTVertex,VertexBase>
describeThePEGAbstractFFSTVertex("ThePEG::AbstractFFSTVertex", "libThePEG.so");

void AbstractFFSTVertex::Init() {

  static ClassDocumentation<AbstractFFSTVertex> documentation
    ("The AbstractFFSTVertex class is the base class for all "
     "fermion-fermion-scalar-tensor interactions in ThePEG.");

}

