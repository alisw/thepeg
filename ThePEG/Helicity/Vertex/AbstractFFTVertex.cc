// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractFFTVertex class.
//

#include "AbstractFFTVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;
using namespace Helicity;

AbstractNoPIOClassDescription<AbstractFFTVertex> 
AbstractFFTVertex::initAbstractFFTVertex;
// Definition of the static class description member.

void AbstractFFTVertex::Init() {

  static ClassDocumentation<AbstractFFTVertex> documentation
    ("The AbstractFFTVertex class is the base class for all fermion-fermion-tensor"
     " interactions in ThePEG");

}

