// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractSSSSVertex class.
//

#include "AbstractSSSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<AbstractSSSSVertex> 
AbstractSSSSVertex::initAbstractSSSSVertex;
// Definition of the static class description member.

void AbstractSSSSVertex::Init() {

  static ClassDocumentation<AbstractSSSSVertex> documentation
    ("The AbstractSSSSVertex class is the base class for all "
     "scalar-scalar-scalar-scalar interactions");

}

