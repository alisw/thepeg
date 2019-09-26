// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVVVVertex class.
//

#include "AbstractVVVVVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractVVVVVertex,VertexBase>
describeThePEGAbstractVVVVVertex("ThePEG::AbstractVVVVVertex", "libThePEG.so");

void AbstractVVVVVertex::Init() {

  static ClassDocumentation<AbstractVVVVVertex> documentation
    ("The AbstractVVVVVertex class is the base class for "
     "vector-vector-vector-vector interactions in ThePEG.");

}

