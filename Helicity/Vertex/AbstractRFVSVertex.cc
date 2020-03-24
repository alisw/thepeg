// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractRFVSVertex class.
//

#include "AbstractRFVSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace ThePEG;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractRFVSVertex,VertexBase>
describeThePEGAbstractRFVSVertex("ThePEG::AbstractRFVSVertex", "libThePEG.so");

void AbstractRFVSVertex::Init() {

  static ClassDocumentation<AbstractRFVSVertex> documentation
    ("The AbstractRFSVertex class is an abstract base class for"
     " the implementation of all spin-3/2 fermion-fermion-vector-scalar vertices");

}

