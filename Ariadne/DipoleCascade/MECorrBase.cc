// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MECorrBase class.
//

#include "MECorrBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MECorrBase.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne;

MECorrBase::~MECorrBase() {}

double MECorrBase::
preweight(tcEmiPtr dipole, tcDipoleStatePtr state, long id,
	  const EmissionType & type) const {
  return 1.0;
}

double MECorrBase::
reweight(tcEmiPtr dip, tcDipoleStatePtr state, long id,
	 Energy2 pt2, vector<double> & genVar,
	 const EmissionType & type) const {
  return 1.0;
}

bool MECorrBase::
finalVeto(tcEmiPtr dipole, tcDipoleStatePtr state, tcParPtr parton,
	  Energy2 pt2, const EmissionType & type) const {
  return false;
}


void MECorrBase::persistentOutput(PersistentOStream & os) const {}

void MECorrBase::persistentInput(PersistentIStream & is, int) {}

AbstractClassDescription<MECorrBase> MECorrBase::initMECorrBase;
// Definition of the static class description member.

void MECorrBase::Init() {

  static ClassDocumentation<MECorrBase> documentation
    ("There is no documentation for the MECorrBase class");

}

