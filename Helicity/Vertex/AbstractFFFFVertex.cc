// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractFFFFVertex class.
//

#include "AbstractFFFFVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractFFFFVertex,VertexBase>
  describeHelicityAbstractFFFFVertex("ThePEG::AbstractFFFFVertex", "libThePEG.so");

void AbstractFFFFVertex::Init() {

  static ClassDocumentation<AbstractFFFFVertex> documentation
    ("There is no documentation for the AbstractFFFFVertex class");

}

