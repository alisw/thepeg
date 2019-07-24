// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractRFVVertex class.
//

#include "AbstractRFVVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractRFVVertex,VertexBase>
describeThePEGAbstractRFVVertex("ThePEG::AbstractRFVVertex", "libThePEG.so");

void AbstractRFVVertex::Init() {

  static ClassDocumentation<AbstractRFVVertex> documentation
    ("The AbstractRFVVertex class provides a base class for all"
     " spin-3/2 fermion-fermion-vector interactions");

}
