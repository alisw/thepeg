// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RRGluonEmission class.
//

#include "RRGluonEmission.h"

#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

RRGluonEmission::~RRGluonEmission() {}

ClonePtr RRGluonEmission::clone() const {
  return new_ptr(*this);
}

ClonePtr RRGluonEmission::fullclone() const {
  return new_ptr(*this);
}

void RRGluonEmission::persistentOutput(PersistentOStream & os) const {
  os << rem1 << rem3 << ounit(mh2, GeV2) << ounit(oph, GeV) << ounit(ophr, GeV)
     << ounit(ophr1, GeV) << ounit(ophr3, GeV)
     << ounit(opr1, GeV) << ounit(opr3, GeV) << ounit(ph, GeV);
}

void RRGluonEmission::persistentInput(PersistentIStream & is, int) {
  is >> rem1 >> rem3 >> iunit(mh2, GeV2) >> iunit(oph, GeV) >> iunit(ophr, GeV)
     >> iunit(ophr1, GeV) >> iunit(ophr3, GeV)
     >> iunit(opr1, GeV) >> iunit(opr3, GeV) >> iunit(ph, GeV);
}

DescribeClass<RRGluonEmission,FSGluonEmission>
describeAriadne5RRGluonEmission("Ariadne5::RRGluonEmission",
				     "libAriadne5.so");

void RRGluonEmission::Init() {}

