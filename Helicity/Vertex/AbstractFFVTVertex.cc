// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractFFVTVertex class.
//

#include "AbstractFFVTVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractFFVTVertex,VertexBase>
describeThePEGAbstractFFVTVertex("ThePEG::AbstractFFVTVertex", "libThePEG.so");

void AbstractFFVTVertex::Init() {

  static ClassDocumentation<AbstractFFVTVertex> documentation
    ("The AbstractFFVTVertex class is the base class for all "
     "fermion-fermion-vector-tensor interactions in ThePEG.");

}

