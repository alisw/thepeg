// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractRFVVertex class.
//

#include "AbstractRFVVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Helicity;

AbstractNoPIOClassDescription<AbstractRFVVertex> 
AbstractRFVVertex::initAbstractRFVVertex;
// Definition of the static class description member.

void AbstractRFVVertex::Init() {

  static ClassDocumentation<AbstractRFVVertex> documentation
    ("The AbstractRFVVertex class provides a base class for all"
     " spin-3/2 fermion-fermion-vector interactions");

}
