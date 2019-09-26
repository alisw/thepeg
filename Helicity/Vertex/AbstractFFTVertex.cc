// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractFFTVertex class.
//

#include "AbstractFFTVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractFFTVertex,VertexBase>
describeThePEGAbstractFFTVertex("ThePEG::AbstractFFTVertex", "libThePEG.so");

void AbstractFFTVertex::Init() {

  static ClassDocumentation<AbstractFFTVertex> documentation
    ("The AbstractFFTVertex class is the base class for all fermion-fermion-tensor"
     " interactions in ThePEG");

}

