// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DYCorr class.
//

#include "DYCorr.h"
#include "ExtendedDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Config/std.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DYCorr.tcc"
#endif


using namespace Ariadne;

DYCorr::~DYCorr() {}

bool DYCorr::canHandle(tcEmiPtr dipole, tcDipoleStatePtr state) const {
  //Make sure there is only one intermediate particle and that is either
  //a Z, W or gamma.
  tPVector intermediates = state->hardSubSys().intermediates();
  if(intermediates.size() != 1){
    return false;
  }
  tPPtr dyboson = intermediates.back();
  if(dyboson->id() != ParticleID::gamma &&
      dyboson->id() != ParticleID::Z0 &&
      dyboson->id() != ParticleID::Wplus &&
      dyboson->id() != ParticleID::Wminus){
    return false;
  }

  //If there are any active particles they should all be children to the
  //Drell-Yan boson.
  HardSubSys::PartonSet active = state->hardSubSys().active();
  tcParticleSet s;
  for(HardSubSys::PartonSet::const_iterator it = active.begin();
      it != active.end(); it++){
    (*it)->getOriginalParents(inserter(s));
  }
  for(tcParticleSet::const_iterator it = s.begin(); it != s.end(); it++){
    tParticleVector p = (*it)->parents();
    for(tParticleVector::const_iterator it2 = p.begin(); it2 != p.end(); it2++){
      if((*it) != dyboson){
        return false;
      }
    }
  }

  //Check that the dipole is an extended dipole between two extended
  //partons with quarks extracted.
  tcExDipPtr edip = dynamic_ptr_cast<tcExDipPtr>(dipole);
  if(!edip || !edip->iSoftRem() || !edip->oSoftRem() ||
      edip->iSoftRem()->isG() || edip->oSoftRem()->isG()){
      return false;
  }

  return true;
}

double DYCorr::reweight(tcEmiPtr dipole, tcDipoleStatePtr state,
    long id, Energy2 pt2, vector<double> & genVar,
    const EmissionType & type) const {
  if(! canHandle(dipole, state)){
    return 1.0;
  }
  if ( Emitter::isType<ExtendedDipole::RRGluon>(type) ){
    double yg = genVar[0];
    Energy2 mh2 = state->hardSubSys().momentum().mass2();
    double yh = state->hardSubSys().momentum().rapidity();
    Energy2 shat = mh2 + 2*pt2 + 2 * sqrt(pt2*(pt2+mh2)) * cosh(yh-yg);
    Energy2 that = -pt2 - sqrt(pt2*(pt2+mh2)) * exp(yh-yg);
    Energy2 uhat = mh2 - shat - that;
    Energy2 S = dynamic_ptr_cast<tcExDipPtr>(dipole)->sdip();
    double x1 = 1.0 - sqrt(pt2/S)*exp( yg);
    double x3 = 1.0 - sqrt(pt2/S)*exp(-yg);
    return pt2/(sqr(x1)+sqr(x3)) * (sqr(that) + sqr(uhat) + 2*shat*mh2) /
      (shat*that*uhat);
  }
  if ( Emitter::isType<ExtendedDipole::ISGtoQQ>(type) ) {
    double xi = genVar[0];
    double z = genVar[1];
    Energy mq = genVar[2]*GeV;
    Energy2 mq2 = sqr(mq);
    Energy mh = state->hardSubSys().momentum().mass();
    Energy2 mh2 = sqr(mh);

    Energy2 shat = mh2 / z;
    Energy2 uhat = -(mq2 * xi + pt2)/(1.0 - xi);
    Energy2 that = mh2 + mq2 - shat - uhat;

    return (sqr(shat) + sqr(uhat) + 2*mh2*that) / sqr(shat) /
      (sqr(z) + sqr(1-z));
  }
  return 1.0;
}

NoPIOClassDescription<DYCorr> DYCorr::initDYCorr;
// Definition of the static class description member.

void DYCorr::Init() {

  static ClassDocumentation<DYCorr> documentation
    ("There is no documentation for the DYCorr class");

}

