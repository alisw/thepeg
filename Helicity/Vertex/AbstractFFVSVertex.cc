// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractFFVSVertex class.
//

#include "AbstractFFVSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace ThePEG;


// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractFFVSVertex,VertexBase>
describeThePEGAbstractFFVSVertex("ThePEG::AbstractFFVSVertex", "libThePEG.so");

void AbstractFFVSVertex::Init() {

  static ClassDocumentation<AbstractFFVSVertex> documentation
    ("The AbstractFFVSVertex class provides a base class for all"
     " fermion-fermionvector-scalar interactions.");

}

