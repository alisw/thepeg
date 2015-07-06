// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourResonanceModel class.
//

#include "ColourResonanceModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

ColourResonanceModel::ColourResonanceModel() {}

ColourResonanceModel::~ColourResonanceModel() {}

IBPtr ColourResonanceModel::clone() const {
  return new_ptr(*this);
}

IBPtr ColourResonanceModel::fullclone() const {
  return new_ptr(*this);
}

PseudoParton ColourResonanceModel::getPseudoParton(tResParPtr p) const {
  // *** ATTENTION *** implement this.
  return PseudoParton(p->isG(), !p->isG(), PseudoParton::colres,
		      p->momentum(), p);
}

bool ColourResonanceModel::stillSpecial(tResParPtr) const {
  // *** Attention *** implement this.
  return true;
}

EmPtr ColourResonanceModel::
generateInternalGluon(const EmitterBase &, const QCDDipole & dip,
		      tResParPtr) const {
  // *** Attention *** implement this.
  return EmPtr();
}

void ColourResonanceModel::
setMomentum(PseudoParton pp, const Lorentz5Momentum & p) const {
  // *** Attention *** implement this.
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ColourResonanceModel::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void ColourResonanceModel::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// The following static variable is needed for the type description
// system in ThePEG.
DescribeClass<ColourResonanceModel,HandlerBase>
describeColourResonanceModel("Ariadne5::ColourResonanceModel", "libAriadne5.so");

void ColourResonanceModel::Init() {

  static ClassDocumentation<ColourResonanceModel> documentation
    ("ColourResonanceModel is a helper class to be used when emitting from "
     "a coloured resonance.");

}

