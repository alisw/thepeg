// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FSGluonEmission class.
//

#include "DipoleSwing.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/EnumIO.h"

#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

DipoleSwing::~DipoleSwing() {}

ClonePtr DipoleSwing::clone() const {
  return new_ptr(*this);
}

ClonePtr DipoleSwing::fullclone() const {
  return new_ptr(*this);
}

void DipoleSwing::setup(QCDDipole & d1, QCDDipole & d2, Time dt) {
  dipoles = make_pair(&d1, &d2);
  radiators.clear();
  radiators.push_back(colourParent = d1.iPart());
  radiators.push_back(antiColourParent = d1.oPart());
  radiators.push_back(d2.iPart());
  radiators.push_back(d2.oPart());
  //    partons = radiators;
  //    swap(partons[1], partons[3]);
  geno = d1.state()->nEmissions() + 1;
  index = d1.colourIndex();
  rho = hbarc/dt;
}

void DipoleSwing::persistentOutput(PersistentOStream & os) const {
  os << dipoles << oenum(forced);
}

void DipoleSwing::persistentInput(PersistentIStream & is, int) {
  is >> dipoles >> ienum(forced);
}

DescribeClass<DipoleSwing,Emission>
describeAriadne5DipoleSwing("Ariadne5::DipoleSwing", "libAriadne5.so");

void DipoleSwing::Init() {}

void DipoleSwing::debugme() const {
  Emission::debugme();
  cerr << "Swinging dipoles (colour index " << index << "):" << endl
       << cdipole->state()->index(dipoles.first)
       << " (partons " << cdipole->state()->index(dipoles.first->iPart())
       << " and " << cdipole->state()->index(dipoles.first->oPart())
       << " colour index " << dipoles.first->colourIndex() << ")" << endl
       << cdipole->state()->index(dipoles.second)
       << " (partons " << cdipole->state()->index(dipoles.second->iPart())
       << " and " << cdipole->state()->index(dipoles.second->oPart())
       << " colour index " << dipoles.second->colourIndex() << ")" << endl;
}
