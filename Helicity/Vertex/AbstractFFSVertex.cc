// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractFFSVertex class.
//

#include "AbstractFFSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Helicity;

AbstractNoPIOClassDescription<AbstractFFSVertex> 
AbstractFFSVertex::initAbstractFFSVertex;
// Definition of the static class description member.

void AbstractFFSVertex::Init() {

  static ClassDocumentation<AbstractFFSVertex> documentation
    ("The AbstractFFSVertex class is an abstract base class for"
     " the implementation of all fermion-fermion-scalar vertices");

}

