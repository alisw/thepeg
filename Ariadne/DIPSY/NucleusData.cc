// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NucleusData class.
//

#include "NucleusData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

NucleusData::NucleusData() : theA(2), theZ(1) {}

NucleusData::~NucleusData() {}

NucleusData::NucleusData(long newId, string newPDGName)
  : ParticleData(newId, newPDGName), theA(2), theZ(1) {}

PDPtr NucleusData::
Create(long newId, string newPDGName) {
  return new_ptr(NucleusData(newId, newPDGName));
}

PDPair NucleusData::
Create(long newId, string newPDGName, string newAntiPDGName) {
  PDPair pap;
  pap.first = new_ptr(NucleusData(newId, newPDGName));
  pap.second = new_ptr(NucleusData(-newId, newAntiPDGName));
  antiSetup(pap);
  return pap;
}

void NucleusData::readSetup(istream & is) {
  ParticleData::readSetup(is);
  is >> theA >> theZ;
}

PDPtr NucleusData::pdclone() const {
  return new_ptr(*this);
}




IBPtr NucleusData::clone() const {
  return new_ptr(*this);
}

IBPtr NucleusData::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NucleusData::persistentOutput(PersistentOStream & os) const {
  os << theA << theZ;
}

void NucleusData::persistentInput(PersistentIStream & is, int) {
  is >> theA >> theZ;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<NucleusData,ThePEG::ParticleData>
  describeDIPSYNucleusData("DIPSY::NucleusData", "libAriadne5.so libDIPSY.so");

void NucleusData::Init() {

  static ClassDocumentation<NucleusData> documentation
    ("There is no documentation for the NucleusData class");

  static Parameter<NucleusData,unsigned int> interfaceA
    ("A",
     "The atomic number for this nucleus.",
     &NucleusData::theA, 2, 2, 0,
     true, false, Interface::lowerlim);

  static Parameter<NucleusData,int> interfaceZ
    ("Z",
     "The number of protons in this nucleus.",
     &NucleusData::theZ, 1, 0, 0,
     true, false, Interface::nolimits);

}

