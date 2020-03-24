// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractRFSSVertex.h class.
//

#include "AbstractRFSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace ThePEG;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractRFSSVertex,VertexBase>
describeThePEGAbstractRFSSVertex("ThePEG::AbstractRFSSVertex", "libThePEG.so");

void AbstractRFSSVertex::Init() {

  static ClassDocumentation<AbstractRFSSVertex> documentation
    ("The AbstractRFSVertex class is an abstract base class for"
     " the implementation of all spin-3/2 fermion-fermion-vector-scalar-scalar vertices");

}

