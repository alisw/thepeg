// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVVSVertex class.
//

#include "AbstractVVVSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

DescribeAbstractNoPIOClass<AbstractVVVSVertex,VertexBase>
describeAbstractVVVVertex("Helicity::AbstractVVVSVertex","libThePEG.so");

void AbstractVVVSVertex::Init() {

  static ClassDocumentation<AbstractVVVSVertex> documentation
    ("The AbstractVVVSVertex class provides the base class for"
     " all vector-vector-vector-scalar interactions in ThePEG."
     "This type of interaction does not occur in normal renormalisable"
     " theories, however it is common in effective theories in particular"
     " that for the Higgs in the infinite top quark mass limit.");

}

