// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Helicity:AbstractFFSSVertex class.
//

#include "AbstractFFSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace ThePEG;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<AbstractFFSSVertex,VertexBase>
describeThePEGAbstractFFSSVertex("ThePEG::AbstractFFSSVertex",
				 "libThePEG.so");

void AbstractFFSSVertex::Init() {

  static ClassDocumentation<AbstractFFSSVertex> documentation
    ("The AbstractFFSSVertex class implements the interaction between a"
     " fermikon-antifermion and two scalars.");

}

