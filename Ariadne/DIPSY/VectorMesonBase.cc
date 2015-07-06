// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonBase class.
//

#include "VectorMesonBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/ParticleData.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

VectorMesonBase::~VectorMesonBase() {}

Energy2 VectorMesonBase::psi2(InvEnergy r, double z) {
  Energy2 sum = 0.0*GeV2;
  for ( int f = 1; f <= abs(maxFlav()); ++f )
    sum += sqr(psi(0, 1, 1, f, r, z)) + sqr(psi(1, 1, -1, f, r, z)) +
      sqr(psi(1, -1, 1, f, r, z)) + sqr(psi(1, 1, 1, f, r, z));
  return sum;
}

void VectorMesonBase::doinit() throw(InitException) {
  WaveFunction::doinit();
  for ( int f = 1; f <= abs(maxFlav()); ++f )
    if ( qmass[f] < 0.0*GeV ) qmass[f] = generator()->getParticleData(f)->mass();
}

void VectorMesonBase::persistentOutput(PersistentOStream & os) const {
  os << theMaxFlav << ounit(qmass, GeV);
}

void VectorMesonBase::persistentInput(PersistentIStream & is, int) {
  is >> theMaxFlav >> iunit(qmass, GeV);
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeAbstractClass<VectorMesonBase,DIPSY::WaveFunction>
  describeDIPSYVectorMesonBase("DIPSY::VectorMesonBase", "libAriadne5.so libDIPSY.so");


void VectorMesonBase::Init() {

  static ClassDocumentation<VectorMesonBase> documentation
    ("The VectorMesonBase class inherits from WaveFunction and is the "
     "base class of all vector meson wave functions. It includes abstract "
     "functions for different polarization and helicity components.");

  static Parameter<VectorMesonBase,int> interfaceMaxFlav
    ("MaxFlav",
     "The maxumim number of flavours considered. If negative only the "
     "corresponding flavour will be considered.",
     &VectorMesonBase::theMaxFlav, 4, -6, 6,
     true, false, Interface::limited);

  static ParVector<VectorMesonBase,Energy> interfaceQuarkMasses
    ("QuarkMasses",
     "The quark masses to be used (zero'th component is always ignored). "
     "If negative, the value in the default corresponding ParticleData "
     "object is used instead.",
     &VectorMesonBase::qmass, GeV, 7, -1.0*GeV, 0*GeV, 0*GeV,
     true, false, Interface::nolimits);

}

