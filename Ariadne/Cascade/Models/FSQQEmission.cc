// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FSQQEmission class.
//

#include "FSQQEmission.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

FSQQEmission::~FSQQEmission() {}

ClonePtr FSQQEmission::clone() const {
  return new_ptr(*this);
}

ClonePtr FSQQEmission::fullclone() const {
  return new_ptr(*this);
}

void FSQQEmission::persistentOutput(PersistentOStream & os) const {
  os << x1 << x3 << ounit(mq, GeV) << ifl << yo << od;
}

void FSQQEmission::persistentInput(PersistentIStream & is, int) {
  is >> x1 >> x3 >> iunit(mq, GeV) >> ifl >> yo >> od;
}

DescribeClass<FSQQEmission,Emission>
describeAriadne5FSQQEmission("Ariadne5::FSQQEmission", "libAriadne5.so");

void FSQQEmission::Init() {}

