// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RRGluonEmitter class.
//

#include "FSGluonEmission.h"
#include "RRGluonEmitter.h"
#include "FSGluonEmitter.h"
#include "SpecialGluonEmitter.h"
#include "RemnantGluonEmitter.h"
#include "RemnantModel.h"
#include "ISGtoQEmitter.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "Ariadne/Cascade/AriadneHandler.h"
#include "Ariadne/Cascade/Resonance.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/DebugItem.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Ariadne/Config/UnitFO.h"


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

RRGluonEmitter::RRGluonEmitter() {}

RRGluonEmitter::~RRGluonEmitter() {}

IBPtr RRGluonEmitter::clone() const {
  return new_ptr(*this);
}

IBPtr RRGluonEmitter::fullclone() const {
  return new_ptr(*this);
}


bool RRGluonEmitter::canHandle(const DipoleBase & e) const {
  tcQCDPtr d = dynamic_cast<const QCDDipole *>(&e);
  if ( !d ) return false;
  tRemParPtr rem1 = dynamic_ptr_cast<tRemParPtr>(d->iPart());
  tRemParPtr rem3 = dynamic_ptr_cast<tRemParPtr>(d->oPart());
  if ( !rem1 || rem1->hard() || !rem3 || rem3->hard() ) return false;
  return true;
}

bool RRGluonEmitter::overrides(const EmitterBase & em, DipoleBase &) const {
  if ( typeid(em) == typeid(FSGluonEmitter) ) return true;
  if ( typeid(em) == typeid(SpecialGluonEmitter) ) return true;
  if ( typeid(em) == typeid(RemnantGluonEmitter) ) return true;
  return false;
}

EmPtr RRGluonEmitter::
generate(const DipoleBase & dipole, Energy rhomin, Energy rhomax) const {
  const QCDDipole & d = dynamic_cast<const QCDDipole &>(dipole);
  tRemParPtr rem1 = dynamic_ptr_cast<tRemParPtr>(d.iPart());
  tRemParPtr rem3 = dynamic_ptr_cast<tRemParPtr>(d.oPart());

  double C = (d.iPart()->isG() || d.oPart()->isG()? 3.0/4.0: 2.0/3.0)/
    Constants::pi;
  int n1 =  d.iPart()->isG()? 3: 2;
  int n3 =  d.oPart()->isG()? 3: 2;
  Energy2 S = d.sdip();
  if ( S <= sqr(2.0*rhomin) ) return EmPtr();
  Energy W = sqrt(S);
  RRGluonEmission e(*this, dipole, rem1, rem3, d.state()->hadronicMomentum());

  C *= preweight(e);

  rhomax = min(rhomax, W/2.0);
  double yint = 2.0*acosh(0.5*W/rhomin);

  while (true) {
    if ( rhomax <= rhomin ) return EmPtr();
    double weight = gluonEngineRR(e, rhomin, rhomax, W, C, yint, n1, n3);
    if ( weight < 0.0 ) return EmPtr();
    rhomax = e.rho;
    if ( weight == 0.0 ) continue;
    weight *= reweight(e);
    if ( weight > UseRandom::rnd() ) break;
  }

  e.colourParent = d.iPart();
  e.antiColourParent = d.oPart();
  e.ymax = log(W/Current<AriadneHandler>()->pTCut());
  e.radiators.push_back(rem1);
  e.radiators.push_back(rem3);

  return new_ptr(e);

}

bool RRGluonEmitter::
perform(const Emission & emission) const {
  return performRRG(dynamic_cast<const RRGluonEmission &>(emission));
}

void RRGluonEmitter::revert(const Emission & emission) const {
  revertRRG(dynamic_cast<const RRGluonEmission &>(emission));
}


double RRGluonEmitter::
gluonEngineRR(RRGluonEmission & e, Energy rhomin, Energy rhomax, Energy W,
	    double C, double yint, int n1, int n3) {
  double weight =
    FSGluonEmitter::gluonEngine(e, rhomin, rhomax, W, C, yint, n1, n3);

  if ( weight < 0.0 ) return weight;

  e.mainParent = e.mainParent = e.rem1;
  if ( sqr(e.x1) > UseRandom::rnd()*(sqr(e.x1) + sqr(e.x3)) )
    e.mainParent = e.rem3;
  try {

    // First generate momenta as if this was a normal emission from
    // the remnants.      
    double phi = 2.0*Constants::pi*UseRandom::rnd();
    e.genmom = getMomenta(sqr(W), e.x1, e.x3, ZERO, ZERO, ZERO, true, true, phi);

    // Now check if the gluon is within the light-cones of the
    // Drell-Yan particle. In that case all recoil goes to the
    // Drell-Yan particle. Otherwise the recoil will only be taken
    // from one of the remnants.
    const Lorentz5Momentum & pg = e.genmom.second;
    Transverse<Energy> pth(e.ophr.x() - pg.x(), e.ophr.y() - pg.y());
    Energy mt = sqrt(pth.pt2() + e.mh2);
    double expy = sqrt(e.ophr.plus()/e.ophr.minus());
    double yr1 = rapidity(e.genmom.first);
    double yg = rapidity(pg);
    double yr3 = rapidity(e.genmom.third);
    double yh = log(expy);

    if ( yr1 < yh || yh < yr3 ) return 0.0;
    double a1 = max(min((yr1 - yg)/(yr1 - yh), 1.0), 0.0);
    double a3 = max(min((yg - yr3)/(yh - yr3), 1.0), 0.0);
    if ( pg.plus() < mt*expy && pg.minus() < mt/expy ) a1 = a3 = 1.0;
    pth = Transverse<Energy>
      (e.ophr.x() + a1*e.genmom.first.x() + a3*e.genmom.third.x(),
       e.ophr.y() + a1*e.genmom.first.y() + a3*e.genmom.third.y());
    mt = sqrt(pth.pt2() + e.mh2);
    Transverse<Energy> pT1((1.0 - a1)*e.genmom.first.x(),
			   (1.0 - a1)*e.genmom.first.y());
    Transverse<Energy> pT3((1.0 - a3)*e.genmom.third.x(),
			   (1.0 - a3)*e.genmom.third.y());
    Energy Pp =
      e.ophr.plus() + e.genmom.first.plus() + e.genmom.third.plus();
    Energy Pm =
      e.ophr.minus() + e.genmom.first.minus() + e.genmom.third.minus();
    if ( a1 == 1.0 ) {
      Pm -= mt/expy;
      Pp -= mt*expy + pT3.pt2()/Pm;
    } else {
      Pp -= mt*expy;
      Pm -= mt/expy + pT1.pt2()/Pp;
    }
    if ( Pp < ZERO || Pm < ZERO ) return 0.0;

    e.genmom.first = e.Rircm*lightCone5(Pp, pT1.pt2()/Pp, pT1);
    e.genmom.second.transform(e.Rircm);
    e.genmom.third = e.Rircm*lightCone5(pT3.pt2()/Pm, Pm, pT3);
    e.ph = e.Rircm*lightCone5(mt*expy, mt/expy, pth);

    LorentzMomentum pem = e.rem1->getBoost()*e.genmom.second;
    double remw1 =
      e.rem1->softSuppression(e.rho, pem.plus()/(e.opr1 + e.ophr1).plus());
    pem = e.rem3->getBoost()*e.genmom.second;
    double remw3 =
      e.rem3->softSuppression(e.rho, pem.plus()/(e.opr3 + e.ophr3).plus());
    double remw = min(remw1, remw3);
    // If the remnant has acquired a transverse momentum an extra
    // suppression is needed for the subsequent emission of a recoil
    // gluon.  The returned weight is the ratio of the new and the old
    // suppression But also if the transverse momentum of the recoil
    // gluon is above the scale of the emission we must veto.
    if ( (e.rem1->getBoost()*e.genmom.first).perp() >= e.rho ) return 0.0;
    if ( (e.rem3->getBoost()*e.genmom.third).perp() >= e.rho ) return 0.0;
    double recw = e.rem1->recoilWeight(e.genmom.second + e.ph, e.genmom.first);
    if ( recw <= 0.0 ) return 0.0;
    recw *= e.rem3->recoilWeight(e.genmom.second + e.ph, e.genmom.third);
    if ( recw <= 0.0 ) return 0.0;
    recw /= e.rem1->recoilWeight()*e.rem3->recoilWeight();
    if ( recw < 1.0 ) remw *= recw;
    weight *= remw;

  } catch ( ImpossibleKinematics ) {
    return 0.0;
  }
    
  return weight;
    
}


bool RRGluonEmitter::performRRG(const RRGluonEmission & e) {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*e.dipole);
  tRemParPtr rem1 = dynamic_ptr_cast<tRemParPtr>(e.radiators[0]);
  tRemParPtr rem3 = dynamic_ptr_cast<tRemParPtr>(e.radiators[1]);
  if ( rem1 == e.mainParent )
    e.Rh = rem1->getHardTransform(e.oph, e.ph);
  else
    e.Rh = rem3->getHardTransform(e.oph, e.ph);
  d.state()->transformHadronicState(e.Rh);
  tParPtr g = FSGluonEmitter::insertGluon(d, e.mainParent).first;
  rem1->setMomentum(e.genmom.first);
  g->momentum() = e.genmom.second;
  rem3->setMomentum(e.genmom.third);

  g->setVertex(e.mainParent->vertex());

  g->emission(&e);
  e.partons.push_back(g);

  return true;  

}

void RRGluonEmitter::revertRRG(const RRGluonEmission & e) {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*e.dipole);
  tRemParPtr rem1 = dynamic_ptr_cast<tRemParPtr>(e.radiators[0]);
  tRemParPtr rem3 = dynamic_ptr_cast<tRemParPtr>(e.radiators[1]);
  FSGluonEmitter::removeGluon(d, d.oPart(), e.mainParent);
  rem1->setMomentum(e.pold.first);
  rem3->setMomentum(e.pold.second);
  d.state()->transformHadronicState(e.Rh.inverse());
  d.state()->untouchHadronicState();
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void RRGluonEmitter::persistentOutput(PersistentOStream &) const {}

void RRGluonEmitter::persistentInput(PersistentIStream &, int) {}

DescribeClass<RRGluonEmitter,ISGtoQEmitter>
describeAriadne5RRGluonEmitter("Ariadne5::RRGluonEmitter", "libAriadne5.so");

void RRGluonEmitter::Init() {

  static ClassDocumentation<RRGluonEmitter> documentation
    ("The RRGluonEmitter class implements the emission of gluons from "
     "a remnant-remnant dipole.");

}

