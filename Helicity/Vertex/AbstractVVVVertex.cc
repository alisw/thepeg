// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractVVVVertex class.
//

#include "AbstractVVVVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

AbstractNoPIOClassDescription<AbstractVVVVertex> 
AbstractVVVVertex::initAbstractVVVVertex;
// Definition of the static class description member.

void AbstractVVVVertex::Init() {

  static ClassDocumentation<AbstractVVVVertex> documentation
    ("The AbstractVVVVertex class provides the base class for"
     " all vector-vector-vector interactions in ThePEG.");

}

