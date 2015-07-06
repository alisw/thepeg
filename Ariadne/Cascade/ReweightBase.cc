// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ReweightBase class.
//

#include "ReweightBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Ariadne5;

Ariadne5::ReweightBase::~ReweightBase() {}

double Ariadne5::ReweightBase::preweight(const Emission &) const {
  return 1.0;
}

double Ariadne5::ReweightBase::reweight(const Emission &) const {
  return 1.0;
}

bool Ariadne5::ReweightBase::finalVeto(const Emission &) const {
  return false;
}

void Ariadne5::ReweightBase::persistentOutput(PersistentOStream & os) const {
  os << mayVeto;
}

void Ariadne5::ReweightBase::persistentInput(PersistentIStream & is, int) {
  is >> mayVeto;
}

DescribeAbstractClass<Ariadne5::ReweightBase,HandlerBase>
describeAriadne5ReweightBase("Ariadne5::ReweightBase", "libAriadne5.so");

void Ariadne5::ReweightBase::Init() {

  static ClassDocumentation<ReweightBase> documentation
    ("ReweightBase is the base class for implementing different ways of "
     "reweighting the dipole emissions in Ariadne.Object of a "
     "concrete sub-class can inserted to a list in "
     "Ariadne5::AriadneHandler and will then be applied to all dipole "
     "emissions.");

}

