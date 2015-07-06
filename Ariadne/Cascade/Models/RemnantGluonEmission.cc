// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RemnantGluonEmission class.
//

#include "RemnantGluonEmission.h"

#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

RemnantGluonEmission::~RemnantGluonEmission() {}

ClonePtr RemnantGluonEmission::clone() const {
  return new_ptr(*this);
}

ClonePtr RemnantGluonEmission::fullclone() const {
  return new_ptr(*this);
}

void RemnantGluonEmission::persistentOutput(PersistentOStream & os) const {
  os << isRecoil;
}

void RemnantGluonEmission::persistentInput(PersistentIStream & is, int) {
  is >> isRecoil;
}

DescribeClass<RemnantGluonEmission,FSGluonEmission>
describeAriadne5RemnantGluonEmission("Ariadne5::RemnantGluonEmission",
				     "libAriadne5.so");

void RemnantGluonEmission::Init() {}

