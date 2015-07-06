// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ISQtoGEmission class.
//

#include "ISQtoGEmission.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

ISQtoGEmission::~ISQtoGEmission() {}

ClonePtr ISQtoGEmission::clone() const {
  return new_ptr(*this);
}

ClonePtr ISQtoGEmission::fullclone() const {
  return new_ptr(*this);
}

void ISQtoGEmission::persistentOutput(PersistentOStream & os) const {
  os << od;
}

void ISQtoGEmission::persistentInput(PersistentIStream & is, int) {
  is >> od;
}

DescribeClass<ISQtoGEmission,ISQEmission>
describeAriadne5ISQtoGEmission("Ariadne5::ISQtoGEmission", "libAriadne5.so");

void ISQtoGEmission::Init() {}

