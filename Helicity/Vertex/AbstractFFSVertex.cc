// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractFFSVertex class.
//

#include "AbstractFFSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractFFSVertex,VertexBase>
describeThePEGAbstractFFSVertex("ThePEG::AbstractFFSVertex", "libThePEG.so");

void AbstractFFSVertex::Init() {

  static ClassDocumentation<AbstractFFSVertex> documentation
    ("The AbstractFFSVertex class is an abstract base class for"
     " the implementation of all fermion-fermion-scalar vertices");

}

