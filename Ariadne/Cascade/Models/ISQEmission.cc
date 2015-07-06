// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ISQEmission class.
//

#include "ISQEmission.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

ISQEmission::~ISQEmission() {}

ClonePtr ISQEmission::clone() const {
  return new_ptr(*this);
}

ClonePtr ISQEmission::fullclone() const {
  return new_ptr(*this);
}

void ISQEmission::persistentOutput(PersistentOStream & os) const {
  os << q << ounit(mq, GeV) << x << exorig << z << xi;
}

void ISQEmission::persistentInput(PersistentIStream & is, int) {
  is >> q >> iunit(mq, GeV) >> x >> exorig >> z >> xi;
}

DescribeClass<ISQEmission,Emission>
describeAriadne5ISQEmission("Ariadne5::ISQEmission", "libAriadne5.so");

void ISQEmission::Init() {}

