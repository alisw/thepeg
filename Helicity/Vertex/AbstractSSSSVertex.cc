// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractSSSSVertex class.
//

#include "AbstractSSSSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractSSSSVertex,VertexBase>
describeThePEGAbstractSSSSVertex("ThePEG::AbstractSSSSVertex", "libThePEG.so");

void AbstractSSSSVertex::Init() {

  static ClassDocumentation<AbstractSSSSVertex> documentation
    ("The AbstractSSSSVertex class is the base class for all "
     "scalar-scalar-scalar-scalar interactions");

}

ScalarWaveFunction AbstractSSSSVertex::evaluate(Energy2,int, tcPDPtr, 
						const ScalarWaveFunction & ,
						const ScalarWaveFunction & ,
						const ScalarWaveFunction & ,
						complex<Energy> ,
						complex<Energy>) {
  assert(false);
}
