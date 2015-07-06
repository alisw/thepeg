// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EMDipole class.
//

#include "EMDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Parton.h"
#include "CascadeHandler.h"
#include "DipoleState.h"
#include "String.h"
#include "MECorrBase.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EMDipole.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Ariadne;

ClonePtr EMDipole::clone() const {
  return new_ptr(*this);
}

bool EMDipole::init(tParPtr ip, tParPtr op, bool respectScale) {
  if(ip->data().charged() and op->data().charged()){
    iPart(ip);
    oPart(op);
    touch();
    lastPT2(-1.0*GeV2);
    setup(sdip()/4.0);
    if ( !respectScale ) return true;
    Energy2 maxscale = sdip()/4.0;
    if ( ip->orig() ) maxscale = min(maxscale, ip->orig()->scale());
    if ( op->orig() ) maxscale = min(maxscale, op->orig()->scale());
    maxScale(maxscale);
    return true;
  }
  return false;
}

Energy2 EMDipole::generate(Energy2 pt2min, Energy2 pt2max) {
  if ( !iPart()->touched() && !oPart()->touched() && !touched() ) {
    return lastPT2();
  }

  lastPT2(-1.0*GeV2);
  double C = preweight(this, state(), ParticleID::gamma, EMPhoton())*0.5/Constants::pi;
  Energy2 S = sdip();
#ifdef Physical_Qty_H
  double Q1 = iPart()->data().charge()/Units::eplus;
  double Q3 = oPart()->data().charge()/Units::eplus;
#else
  double Q1 = iPart()->data().charge();
  double Q3 = oPart()->data().charge();
#endif
  double y1 = max(iPart()->momentum().mass2()/S, 0.0);
  double y3 = max(oPart()->momentum().mass2()/S, 0.0);
  pt2max = min(pt2max, S/4.0);
  if ( pt2max <= pt2min ) return lastPT2();
  double yint = 2.0*acosh(sqrt(S/pt2min)/2.0);

  //Vector containing x1 and x3
  vector<double> var(2);

  while (true) {
    pt2max = rndsudEM(2.0*C*yint*(Q1*Q1 + Q3*Q3), pt2max, pt2min);
    if ( pt2max <= pt2min ) break;
    double ymax = acosh(sqrt(S/pt2max)/2.0);
    double y = 2.0*ymax*handler()->rnd() - ymax;
    var[0] = 1.0 - sqrt(pt2max/S)*exp(-y) + y1 - y3;
    var[1] = 1.0 - sqrt(pt2max/S)*exp( y) + y3 - y1;

    if ( !check(var[0], var[1], y1, 0.0, y3) ) continue;

    double weight = (sqr(var[0]) + sqr(var[1]))/2.0;
    weight *= 2.0*ymax/yint;

    weight *= sqr((Q1*(1 - var[0]) - Q3*(1 - var[1])) /
        (2 - var[0] - var[1])) / (sqr(Q1) + sqr(Q3));

    weight *= reweight(this, state(), ParticleID::gamma, pt2max, var, EMPhoton());

    if ( weight < handler()->rnd() ) continue;

    lastPT2(pt2max);
    lastEType(EMPhoton());
    genVar.swap(var);
    break;
  }
  return lastPT2();
}

tParPtr EMDipole::perform() {
  Lorentz5Momentum p1(iPart()->momentum().mass());
  Lorentz5Momentum p2;
  Lorentz5Momentum p3(oPart()->momentum().mass());
  try {
    SimplePhaseSpace::CMS(p1, p2, p3, sdip(), genVar[0], genVar[1], 0.0, 0.0, 0.0);
  }
  catch ( ImpossibleKinematics ) {
    return tParPtr();
  }
  double Psi = Constants::pi - p3.theta();
  double beta = 0.0;
  if (sqr(genVar[1]) > handler()->rnd()*(sqr(genVar[0]) + sqr(genVar[1])) )
    beta = Psi;

  LorentzRotation R;
  R.rotateY(-beta);
  R.rotateZ(handler()->rnd(2.0*Constants::pi));
  R.transform(Utilities::getBoostFromCM(make_pair(iPart()->momentum(),
						  oPart()->momentum())));
  p1.transform(R);
  p2.transform(R);
  p3.transform(R);


  ParPtr gamma =
    state()->create(handler()->generator()->getParticleData(ParticleID::gamma));
  if ( !gamma ) return gamma;
  gamma->momentum() = p2;
  gamma->parents(make_pair(iPart(), oPart()));

  iPart()->momentum() = p1;
  oPart()->momentum() = p3;

  iPart()->touch();
  oPart()->touch();
  touch();

  return gamma;
}

Energy2 EMDipole::sdip() const {
  return (iPart()->momentum() + oPart()->momentum()).m2();
}

void EMDipole::fillReferences(CloneSet & cset) const {
  CascadeBase::fillReferences(cset);
  cset.insert(oPart());
  cset.insert(iPart());
}

void EMDipole::rebind(const TranslationMap & trans) {
  CascadeBase::rebind(trans);
  oPart(trans.translate(oPart()));
  iPart(trans.translate(iPart()));
}

void EMDipole::persistentOutput(PersistentOStream & os) const {
  os << theIPart << theOPart;
}

void EMDipole::persistentInput(PersistentIStream & is, int) {
  is >> theIPart >> theOPart;
}

ClassDescription<EMDipole> EMDipole::initEMDipole;
// Definition of the static class description member.

void EMDipole::Init() {

  static ClassDocumentation<EMDipole> documentation
    ("There is no documentation for the EMDipole class");

}

void EMDipole::debugme() const {
  Emitter::debugme();
}

