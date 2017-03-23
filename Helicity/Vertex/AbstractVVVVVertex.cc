// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVVVVertex class.
//

#include "AbstractVVVVVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<AbstractVVVVVertex> 
AbstractVVVVVertex::initAbstractVVVVVertex;
// Definition of the static class description member.

void AbstractVVVVVertex::Init() {

  static ClassDocumentation<AbstractVVVVVertex> documentation
    ("The AbstractVVVVVertex class is the base class for "
     "vector-vector-vector-vector interactions in ThePEG.");

}

