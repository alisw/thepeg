// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVSSVertex class.
//

#include "AbstractVVSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<AbstractVVSSVertex> 
AbstractVVSSVertex::initAbstractVVSSVertex;
// Definition of the static class description member.

void AbstractVVSSVertex::Init() {

  static ClassDocumentation<AbstractVVSSVertex> documentation
    ("The AbstractVVSSVertex class is the base class for "
     "vector-vector-scalar-scalar interactions in ThePEG");

}

