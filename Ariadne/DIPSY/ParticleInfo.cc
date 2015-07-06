// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ParticleInfo class.
//

#include "ParticleInfo.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

ParticleInfo::ParticleInfo(tcPartonPtr p): theParton(p) {}

ParticleInfo::~ParticleInfo() {}

tcPartonPtr ParticleInfo::getParton(const Particle & particle) {
  for ( int i = 0, N = particle.getInfo().size(); i < N; ++i )
    if ( const ParticleInfo * pin =
	 dynamic_cast<const ParticleInfo *>(particle.getInfo()[i].operator->()) )
      return pin->parton();
  return tcPartonPtr();
}


void ParticleInfo::persistentOutput(PersistentOStream & os) const {
  os << theParton;
}

void ParticleInfo::persistentInput(PersistentIStream & is, int) {
  is >> theParton;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<ParticleInfo,ThePEG::EventInfoBase>
  describeDIPSYParticleInfo("DIPSY::ParticleInfo", "ParticleInfo.so");

void ParticleInfo::Init() {

  static ClassDocumentation<ParticleInfo> documentation
    ("There is no documentation for the ParticleInfo class");

}

