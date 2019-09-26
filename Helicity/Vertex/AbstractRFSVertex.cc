// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractRFSVertex class.
//

#include "AbstractRFSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractRFSVertex,VertexBase>
describeThePEGAbstractRFSVertex("ThePEG::AbstractRFSVertex", "libThePEG.so");

void AbstractRFSVertex::Init() {

  static ClassDocumentation<AbstractRFSVertex> documentation
    ("The AbstractRFSVertex class is an abstract base class for"
     " the implementation of all spin-3/2 fermion-fermion-scalar vertices");

}

