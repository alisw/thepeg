// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleBase class.
//

#include "DipoleBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/ColourLine.h"

using namespace Ariadne5;

ColinePtr DipoleBase::colourLine() const {
  return ColinePtr();
}

bool DipoleBase::checkIntegrity() {
  return true;
}

void DipoleBase::persistentOutput(PersistentOStream & os) const {
  os << ounit(theRhoCut, GeV);
}

void DipoleBase::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theRhoCut, GeV);
}

// Definition of the static class description member.
DescribeAbstractClass<DipoleBase,CascadeBase>
describeAriadne5DipoleBase("Ariadne5::DipoleBase", "libAriadne5.so");

void DipoleBase::Init() {}

