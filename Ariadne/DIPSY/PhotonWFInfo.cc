// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhotonWFInfo class.
//

#include "PhotonWFInfo.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PhotonWFInfo.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

PhotonWFInfo::~PhotonWFInfo() {}

void PhotonWFInfo::persistentOutput(PersistentOStream & os) const {
  os << theZ << thePol << theH << theHbar << theFlav;
}

void PhotonWFInfo::persistentInput(PersistentIStream & is, int) {
  is >> theZ >> thePol >> theH >> theHbar >> theFlav;
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<PhotonWFInfo,DIPSY::WFInfo>
  describeDIPSYPhotonWFInfo("DIPSY::PhotonWFInfo", "libAriadne5.so libDIPSY.so");

void PhotonWFInfo::Init() {}

