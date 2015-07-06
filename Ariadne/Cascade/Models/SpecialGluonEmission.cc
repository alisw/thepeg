// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SpecialGluonEmission class.
//

#include "SpecialGluonEmission.h"
#include "Ariadne/Cascade/Parton.h"
#include "ThePEG/Utilities/EnumIO.h"

#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

SpecialGluonEmission::~SpecialGluonEmission() {}

ClonePtr SpecialGluonEmission::clone() const {
  return new_ptr(*this);
}

ClonePtr SpecialGluonEmission::fullclone() const {
  return new_ptr(*this);
}

void SpecialGluonEmission::persistentOutput(PersistentOStream & os) const {
  os << pip.isGluon << pip.noRecoil << oenum(pip.type)
     << ounit(pip.p, GeV) << pip.realParton
     << pop.isGluon << pop.noRecoil << oenum(pop.type)
     << ounit(pop.p, GeV) << pop.realParton
     << ounit(S, GeV2) << wrem1 << wrem3;
}

void SpecialGluonEmission::persistentInput(PersistentIStream & is, int) {
  is >> pip.isGluon >> pip.noRecoil >> ienum(pip.type)
     >> iunit(pip.p, GeV) >> pip.realParton
     >> pop.isGluon >> pop.noRecoil >> ienum(pop.type)
     >> iunit(pop.p, GeV) >> pop.realParton
     >> iunit(S, GeV2) >> wrem1 >> wrem3;
}

DescribeClass<SpecialGluonEmission,FSGluonEmission>
describeAriadne5SpecialGluonEmission("Ariadne5::SpecialGluonEmission",
				     "libAriadne5.so");

void SpecialGluonEmission::Init() {}

