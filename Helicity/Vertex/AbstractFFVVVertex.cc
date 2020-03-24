// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractFFVVVertex class.
//

#include "AbstractFFVVVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractFFVVVertex,VertexBase>
describeThePEGAbstractFFVVVertex("ThePEG::AbstractFFVVVertex", "libThePEG.so");

void AbstractFFVVVertex::Init() {

  static ClassDocumentation<AbstractFFVVVertex> documentation
    ("The AbstractFFVVVertex class provides a base class for all"
     " fermion-fermion-vector-vector interactions");

}

