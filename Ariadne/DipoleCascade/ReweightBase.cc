// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ReweightBase class.
//

#include "ReweightBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ReweightBase.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne;

Ariadne::ReweightBase::~ReweightBase() {}

double Ariadne::ReweightBase::
preweight(tcEmiPtr dipole, tcDipoleStatePtr state, long id,
	  const EmissionType & type) const {
  return 1.0;
}

double Ariadne::ReweightBase::
reweight(tcEmiPtr dip, tcDipoleStatePtr state, long id,
	 Energy2 pt2, vector<double> & genVar,
	 const EmissionType & type) const {
  return 1.0;
}

bool Ariadne::ReweightBase::
finalVeto(tcEmiPtr dipole, tcDipoleStatePtr state, tcParPtr parton,
	  Energy2 pt2, const EmissionType & type) const {
  return false;
}


void Ariadne::ReweightBase::persistentOutput(PersistentOStream & os) const {}

void Ariadne::ReweightBase::persistentInput(PersistentIStream & is, int) {}

AbstractClassDescription<Ariadne::ReweightBase> Ariadne::ReweightBase::initReweightBase;
// Definition of the static class description member.

void Ariadne::ReweightBase::Init() {

  static ClassDocumentation<ReweightBase> documentation
    ("ReweightBase is the base class for implementing different ways of "
     "reweighting the dipole emissions in Ariadne.Object of a "
     "concrete sub-class can inserted to a list in "
     "Ariadne::CascadeHandler and will then be applied to all dipole "
     "emissions, even those which are already correced by an MECorrBase "
     "object.");

}

