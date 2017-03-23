// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractSSSVertex class.
//

#include "AbstractSSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<AbstractSSSVertex> 
AbstractSSSVertex::initAbstractSSSVertex;
// Definition of the static class description member.

void AbstractSSSVertex::Init() {

  static ClassDocumentation<AbstractSSSVertex> documentation
    ("The AbstractSSSVertex class is the base class for all "
     "scalar-scalar-scalar interactions");

}

