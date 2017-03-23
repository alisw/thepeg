// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractRFSVertex class.
//

#include "AbstractRFSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Helicity;

AbstractNoPIOClassDescription<AbstractRFSVertex> 
AbstractRFSVertex::initAbstractRFSVertex;
// Definition of the static class description member.

void AbstractRFSVertex::Init() {

  static ClassDocumentation<AbstractRFSVertex> documentation
    ("The AbstractRFSVertex class is an abstract base class for"
     " the implementation of all spin-3/2 fermion-fermion-scalar vertices");

}

