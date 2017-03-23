// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVVTVertex class.
//

#include "AbstractVVVTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<AbstractVVVTVertex> 
AbstractVVVTVertex::initAbstractVVVTVertex;
// Definition of the static class description member.

void AbstractVVVTVertex::Init() {

  static ClassDocumentation<AbstractVVVTVertex> documentation
    ("The AbstractVVVTVertex class is the base class for vector-vector-vector-tensor"
     " interactions in ThePEG");

}

