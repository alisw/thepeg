// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVVVertex class.
//

#include "AbstractVVVVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractVVVVertex,VertexBase>
describeThePEGAbstractVVVVertex("ThePEG::AbstractVVVVertex", "libThePEG.so");

void AbstractVVVVertex::Init() {

  static ClassDocumentation<AbstractVVVVertex> documentation
    ("The AbstractVVVVertex class provides the base class for"
     " all vector-vector-vector interactions in ThePEG.");

}

