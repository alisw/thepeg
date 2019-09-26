// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVSVertex class.
//

#include "AbstractVVSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractVVSVertex,VertexBase>
describeThePEGAbstractVVSVertex("ThePEG::AbstractVVSVertex", "libThePEG.so");

void AbstractVVSVertex::Init() {

  static ClassDocumentation<AbstractVVSVertex> documentation
    ("The AbstractVVSVertex class is the base class for the "
     "vector-vector-scalar interaction");

}

