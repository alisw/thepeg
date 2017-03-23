// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVTVertex class.
//

#include "AbstractVVTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<AbstractVVTVertex> 
AbstractVVTVertex::initAbstractVVTVertex;
// Definition of the static class description member.

void AbstractVVTVertex::Init() {

  static ClassDocumentation<AbstractVVTVertex> documentation
    ("The AbstractVVTVertex class is the base class for all "
     "vector-vector-tensor interactions in ThEPEG");

}

