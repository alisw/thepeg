// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WaveFunction class.
//

#include "WaveFunction.h"
#include "DipoleEventHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

WaveFunction::~WaveFunction() {}

void WaveFunction::fixValence(Step &, tPPtr, const vector<PPtr> &) const {}

void WaveFunction::initialize(const DipoleEventHandler & eh) {
  theEventHandler = &eh;
}

void WaveFunction::persistentOutput(PersistentOStream & os) const {
  os << theEventHandler << theParticle;
}

void WaveFunction::persistentInput(PersistentIStream & is, int) {
  is >> theEventHandler >> theParticle;
}

void WaveFunction::setParticle(PDPtr p) {
  theParticle = p;
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeAbstractClass<WaveFunction,HandlerBase>
  describeDIPSYWaveFunction("DIPSY::WaveFunction", "libAriadne5.so libDIPSY.so");

void WaveFunction::Init() {

  static ClassDocumentation<WaveFunction> documentation
    ("WaveFunction is the base class for wavefunction objects capable of "
     "generating initial DipoleState objects.");

  static Reference<WaveFunction,ParticleData> interfaceParticle
    ("Particle",
     "The corresponding particle.",
     &WaveFunction::theParticle, true, false, true, false, false,
     &WaveFunction::setParticle,
     (Ptr<ParticleData>::pointer(WaveFunction::*)()const)(0),
     (bool(WaveFunction::*)(Ptr<ParticleData>::const_pointer)const)(0));

}

