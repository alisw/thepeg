// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EMDipole class.
//

#include "EMDipole.h"
#include "Parton.h"
#include "AriadneHandler.h"

#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Ariadne5;

EMDipole::EMDipole()
  : theSystem(0) {
  theRhoCut = Current<AriadneHandler>()->pTCutEM();
}

ClonePtr EMDipole::clone() const {
  return new_ptr(*this);
}

Energy2 EMDipole::sdip() const {
  return (iPart()->momentum() + oPart()->momentum()).m2();
}

void EMDipole::fillReferences(CloneSet & cset) const {
  DipoleBase::fillReferences(cset);
  cset.insert(oPart());
  cset.insert(iPart());
}

void EMDipole::rebind(const TranslationMap & trans) {
  DipoleBase::rebind(trans);
  oPart(trans.translate(oPart()));
  iPart(trans.translate(iPart()));
}

void EMDipole::persistentOutput(PersistentOStream & os) const {
  os << theIPart << theOPart << theSystem;
}

void EMDipole::persistentInput(PersistentIStream & is, int) {
  is >> theIPart >> theOPart >> theSystem;
}

DescribeClass<EMDipole,DipoleBase>
describeAriadne5EMDipole("Ariadne5::EMDipole", "libAriadne5.so");

void EMDipole::Init() {}

void EMDipole::debugme() const {
  DipoleBase::debugme();
}

