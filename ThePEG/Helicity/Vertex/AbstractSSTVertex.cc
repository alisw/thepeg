// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractSSTVertex class.
//

#include "AbstractSSTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<AbstractSSTVertex> 
AbstractSSTVertex::initAbstractSSTVertex;
// Definition of the static class description member.

void AbstractSSTVertex::Init() {

  static ClassDocumentation<AbstractSSTVertex> documentation
    ("The AbstractSSTVertex class is the base class for scalar-scalar-tensor"
     "interactions in ThePEG.");

}

