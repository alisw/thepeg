// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RemnantGluonEmitter class.
//

#include "RemnantGluonEmission.h"
#include "RemnantGluonEmitter.h"
#include "RemnantModel.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "Ariadne/Cascade/AriadneHandler.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

RemnantGluonEmitter::RemnantGluonEmitter() {}

RemnantGluonEmitter::~RemnantGluonEmitter() {}

IBPtr RemnantGluonEmitter::clone() const {
  return new_ptr(*this);
}

IBPtr RemnantGluonEmitter::fullclone() const {
  return new_ptr(*this);
}


bool RemnantGluonEmitter::canHandle(const DipoleBase & e) const {
  if ( const QCDDipole * d = dynamic_cast<const QCDDipole *>(&e) ) {
    if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d->iPart()) )
      if ( !r->hard() ) return true;
    if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d->oPart()) )
      if ( !r->hard() ) return true;
  }
  return false;
}

bool RemnantGluonEmitter::
overrides(const EmitterBase & em, DipoleBase &) const {
  return false;
}

EmPtr RemnantGluonEmitter::
generate(const DipoleBase & dipole, Energy rhomin, Energy rhomax) const {
  const QCDDipole & d = dynamic_cast<const QCDDipole &>(dipole);

  EmSel sel(rhomin);

  if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d.iPart()) )
    if ( !r->hard() )
      sel = generateInitialStateGluon(d, r, rhomin, rhomax);
  if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d.oPart()) )
    if ( !r->hard() )
      sel = generateInitialStateGluon(d, r, rhomin, rhomax);

  return sel;
}


bool RemnantGluonEmitter::
perform(const Emission & emission) const {
  const RemnantGluonEmission & e =
    dynamic_cast<const RemnantGluonEmission &>(emission);
  return performInitialStateGluon(e);
}

void RemnantGluonEmitter::revert(const Emission & emission) const {
  const RemnantGluonEmission & e =
    dynamic_cast<const RemnantGluonEmission &>(emission);
  revertInitialStateGluon(e);
  return;
}

EmPtr RemnantGluonEmitter::
generateInitialStateGluon(const QCDDipole & dip, tRemParPtr rem,
			  Energy rhomin, Energy rhomax) const {
  LorentzRotation Rr = rem->getBoost();
  Lorentz5Momentum ph = rem->state()->hadronicMomentum();
  Lorentz5Momentum pr = rem->momentum();
  Lorentz5Momentum phr = Rr*ph;
  Lorentz5Momentum prr = Rr*pr;
  Energy2 S = (ph + pr).m2();
  if ( S <= sqr(rhomin*2.0) ) return EmPtr();
  Energy W = sqrt(S);
  double y1 = max(pr.mass2()/S, 0.0);
  double y3 = max(ph.mass2()/S, 0.0);
  
  double C = 3.0/(4.0*Constants::pi);

  rhomax = min(rhomax, W/2.0);
  if ( rhomax <= rhomin ) return EmPtr();
  double yint = 2.0*acosh(0.5*W/rhomin);

  RemnantGluonEmission e(*this, dip, y1, y3);

  while (true) {
    double weight =
      FSGluonEmitter::gluonEngine(e, rhomin, rhomax, W, C, yint, 3, 3);
    if ( weight < 0.0 ) return EmPtr();
    rhomax = e.rho;
    if ( weight == 0.0 ) continue;

    try {

      // Generate the momenta in the relevant remnant system
      e.genmom =
	FSGluonEmitter::getMomenta(S, e.x1, e.x3, pr.mass(), ZERO, ph.mass(),
				    true, true,
				    UseRandom::rnd(2.0*Constants::pi), pr, ph);
      // Now calculate the soft suppression due to the fraction of
      // positive momentum taken from the incoming particle.
      LorentzMomentum pem = Rr*e.genmom.second;
      double remw = rem->softSuppression(e.rho, pem.plus()/(prr + phr).plus());
      // If the remnant has acquired a transverse momentum an extra
      // suppression is needed for the subsequent emission of a recoil
      // gluon.  The returned weight is the ratio of the new and the old
      // suppression But also if the transverse momentum of the recoil
      // gluon is above the scale of the emission we must veto.
      LorentzMomentum prem = Rr*e.genmom.first;
      if ( prem.perp() >= e.rho ) continue;
      double recw =
	rem->recoilWeight(e.genmom.third + e.genmom.second, e.genmom.first);
      if ( recw <= 0.0 ) continue;
      recw /= rem->recoilWeight();
      if ( recw < 1.0 ) remw *= recw;
      weight *= remw;
 
    } catch ( ImpossibleKinematics ) {
      continue;
    }

    if ( weight > UseRandom::rnd() ) break;

  }

  e.colourParent = dip.iPart();
  e.antiColourParent = dip.oPart();  
  e.mainParent = rem;
  e.pold = make_pair(pr, ph);
  e.ymax = log(W/Current<AriadneHandler>()->pTCut());
  e.radiators.push_back(rem);

  return new_ptr(e);

}

bool RemnantGluonEmitter::
performInitialStateGluon(const RemnantGluonEmission & e) const {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*e.dipole);
  tRemParPtr rem = dynamic_ptr_cast<tRemParPtr>(e.mainParent);
  LorentzMomentum phold = d.state()->hadronicMomentum();
  LorentzMomentum phnew = e.genmom.third;
  e.Rh = rem->getHardTransform(phold, phnew);
  d.state()->transformHadronicState(e.Rh);
  tParPtr g = FSGluonEmitter::insertGluon(d, rem).first;
  g->momentum() = e.genmom.second;
  rem->setMomentum(e.genmom.first);

  g->setVertex(rem->vertex());

  g->emission(&e);
  e.partons.push_back(g);

  return true;
}

void RemnantGluonEmitter::
revertInitialStateGluon(const RemnantGluonEmission & e) const {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*e.dipole);
  tRemParPtr rem = dynamic_ptr_cast<tRemParPtr>(e.mainParent);
  FSGluonEmitter::removeGluon(d, d.oPart(), rem);
  //  LorentzMomentum phnew = e.genmom.third;
  //  LorentzMomentum phold = e.pold.second;
  rem->setMomentum(e.pold.first);
  d.state()->transformHadronicState(e.Rh.inverse());
  d.state()->untouchHadronicState();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void RemnantGluonEmitter::persistentOutput(PersistentOStream &) const {}

void RemnantGluonEmitter::persistentInput(PersistentIStream &, int) {}

DescribeClass<RemnantGluonEmitter,FSGluonEmitter>
describeAriadne5RemnantGluonEmitter("Ariadne5::RemnantGluonEmitter",
				    "libAriadne5.so");

void RemnantGluonEmitter::Init() {

  static ClassDocumentation<RemnantGluonEmitter> documentation
    ("The RemnantGluonEmitter class implements the soft gluon "
     "radiation from a remnant dipole.");

}

