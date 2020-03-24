// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractRFVVVertex class.
//

#include "AbstractRFVVVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace ThePEG;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractRFVVVertex,VertexBase>
describeThePEGAbstractRFVVVertex("ThePEG::AbstractRFVVVertex", "libThePEG.so");

void AbstractRFVVVertex::Init() {

  static ClassDocumentation<AbstractRFVVVertex> documentation
    ("The AbstractRFSVertex class is an abstract base class for"
     " the implementation of all spin-3/2 fermion-fermion-vector-vector vertices");

}

