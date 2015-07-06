// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractFFVTVertex class.
//

#include "AbstractFFVTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<AbstractFFVTVertex> 
AbstractFFVTVertex::initAbstractFFVTVertex;
// Definition of the static class description member.

void AbstractFFVTVertex::Init() {

  static ClassDocumentation<AbstractFFVTVertex> documentation
    ("The AbstractFFVTVertex class is the base class for all "
     "fermion-fermion-vector-tensor interactions in ThePEG.");

}

