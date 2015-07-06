// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVSVertex class.
//

#include "AbstractVVSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<AbstractVVSVertex> 
AbstractVVSVertex::initAbstractVVSVertex;
// Definition of the static class description member.

void AbstractVVSVertex::Init() {

  static ClassDocumentation<AbstractVVSVertex> documentation
    ("The AbstractVVSVertex class is the base class for the "
     "vector-vector-scalar interaction");

}

