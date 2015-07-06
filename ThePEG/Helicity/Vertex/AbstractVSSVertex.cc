// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVSSVertex class.
//

#include "AbstractVSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<AbstractVSSVertex> 
AbstractVSSVertex::initAbstractVSSVertex;
// Definition of the static class description member.

void AbstractVSSVertex::Init() {

  static ClassDocumentation<AbstractVSSVertex> documentation
    ("The AbstractVSSVertex class is the base class for all "
     "vector-scalar-scalar interactions in ThePEG.");

}

