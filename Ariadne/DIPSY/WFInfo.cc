// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WFInfo class.
//

#include "WFInfo.h"
#include "WaveFunction.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

WFInfo::~WFInfo() {}

tcWFInfoPtr WFInfo::getWFInfo(const Particle & particle) {
  for ( int i = 0, N = particle.getInfo().size(); i < N; ++i )
    if ( const WFInfo * win =
	 dynamic_cast<const WFInfo *>(particle.getInfo()[i].operator->()) )
      return win;
  return tcWFInfoPtr();
}

void WFInfo::persistentOutput(PersistentOStream & os) const {
  os << theWF << ounit(theR, InvGeV);
}

void WFInfo::persistentInput(PersistentIStream & is, int) {
  is >> theWF >> iunit(theR, InvGeV);
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<WFInfo,ThePEG::EventInfoBase>
  describeDIPSYWFInfo("DIPSY::WFInfo", "libAriadne5.so libDIPSY.so");

void WFInfo::Init() {}

