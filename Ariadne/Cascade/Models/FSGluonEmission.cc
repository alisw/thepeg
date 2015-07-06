// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FSGluonEmission class.
//

#include "FSGluonEmission.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

FSGluonEmission::~FSGluonEmission() {}

ClonePtr FSGluonEmission::clone() const {
  return new_ptr(*this);
}

ClonePtr FSGluonEmission::fullclone() const {
  return new_ptr(*this);
}

void FSGluonEmission::persistentOutput(PersistentOStream & os) const {
  os << x1 << x3 << y1 << y3;
}

void FSGluonEmission::persistentInput(PersistentIStream & is, int) {
  is >> x1 >> x3 >> y1 >> y3;
}

DescribeClass<FSGluonEmission,Emission>
describeAriadne5FSGluonEmission("Ariadne5::FSGluonEmission", "libAriadne5.so");

void FSGluonEmission::Init() {}

