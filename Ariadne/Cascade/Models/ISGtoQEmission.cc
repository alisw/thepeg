// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ISGtoQEmission class.
//

#include "ISGtoQEmission.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

ISGtoQEmission::~ISGtoQEmission() {}

ClonePtr ISGtoQEmission::clone() const {
  return new_ptr(*this);
}

ClonePtr ISGtoQEmission::fullclone() const {
  return new_ptr(*this);
}

DescribeNoPIOClass<ISGtoQEmission,ISQEmission>
describeAriadne5ISGtoQEmission("Ariadne5::ISGtoQEmission", "libAriadne5.so");

void ISGtoQEmission::Init() {}

