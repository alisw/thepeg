// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RemnantModel class.
//

#include "RemnantModel.h"
#include "FSGluonEmitter.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "Ariadne/Cascade/AriadneHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/UtilityBase.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

RemnantModel::RemnantModel() {}

RemnantModel::~RemnantModel() {}

IBPtr RemnantModel::clone() const {
  return new_ptr(*this);
}

IBPtr RemnantModel::fullclone() const {
  return new_ptr(*this);
}

PseudoParton RemnantModel::getPseudoParton(tRemParPtr rem) const {
  return PseudoParton(rem->isG(), true, PseudoParton::remnant,
		      rem->momentum(), rem);
}

double RemnantModel::
reweightFS(tcRemParPtr rem, Energy rho, LorentzMomentum pem,
	   LorentzMomentum prem, tcPDPtr extracted) const {
  // Boost momenta of the emitted gluon and the hard subsystem so that
  // the (original) remnant is along the +z-axis. Veto if the emitted
  // gluon has more negative momentum than the original hard
  // subsystem.
  LorentzMomentum premold = rem->getBoost()*rem->momentum();
  LorentzMomentum ph = rem->getBoost()*rem->state()->hadronicMomentum();
  pem.transform(rem->getBoost());

  // If the emission is not inside the the momentum region of the hard
  // subsystem, then it is not considered to be a final-state emission
  // and is vetoed. In case we have a hard remnant as the only parton
  // in the hard sybsystem, then we assume DIS and the momentum region
  // is given by Q2 (which is the same as the square of mu()).
  Energy phminus = ph.minus();
  if ( rem->state()->hadronicFS().size() == 1 ) {
    tRemParPtr hrem =
      dynamic_ptr_cast<tRemParPtr>(*(rem->state()->hadronicFS().begin()));
    if ( hrem && hrem->hard() ) 
      phminus = sqr(hrem->mu())/(premold + ph).m();
  }
  if ( pem.minus() > phminus ) return 0.0;

  // Now calculate the soft suppression due to the fraction of
  // positive momentum taken from the incoming particle.
  double weight =
    softSuppression(rho, pem.plus()/(premold + ph).plus(),  rem->x(),
		    rem->mu(), rem->alpha(), rem->beta());

  // If the remnant has acquired a transverse momentum an extra
  // suppression is needed for the subsequent emission of a recoil
  // gluon.  The returned weight is the ratio of the new and the old
  // suppression But also if the transverse momentum of the recoil
  // gluon is above the scale of the emission we must veto.
  prem.transform(rem->getBoost());
  if ( prem.perp() >= rho ) return 0.0;
  double recw = recoilWeight(ph + premold - prem, prem, rem, extracted);
  if ( recw <= 0.0 ) return 0.0;
  recw /= recoilWeight(ph, premold, rem);
  return recw < 1.0? recw*weight: weight;

}

double RemnantModel::
softSuppression(Energy rho, double xplus, double x,
	       Energy mu, double alpha, double spow) const {

  // For a given emission scale rho, the maximum allowed positive
  // momentum to be taken from the incoming particle is:
  if ( xplus < x ) return 1.0;
  double xpmax = pow(mu/rho, alpha);
  if ( xplus <= xpmax ) return 1.0;

  // For a given positive momentum taken from the incoming particle,
  // the maximum allowed emission scale is:
  Energy rhomax = mu*pow(xplus, -1.0/alpha);
  rhomax = max(rhomax, rho*x/xplus);

  // Reweight with the ratio of maximum and actual emission scale to
  // some power.
  return pow(rhomax/rho, spow);

}

double RemnantModel::
recoilWeight(const LorentzMomentum & ph, const LorentzMomentum & pr,
	     tcRemParPtr rem, tcPDPtr extracted) const {
  if ( !extracted ) extracted = &rem->extractedData();
  Energy effMr = effectiveMass(rem, extracted);
  if ( (ph + pr).m2() < sqr(ph.m() + effMr + pr.perp()) ) return 0.0;
  double alpha = rem->alpha();
  double beta = rem->beta();
  Energy mu = rem->mu();
  // Only transverse momenta above the cutoff gives a suppression.
  Energy kt = pr.perp();
  if ( kt <= Current<AriadneHandler>()->pTCut() ) return 1.0;

  // First calculate the limits on the positive momentum of the recoil
  // and return zero if no phase space is available for the recoil
  // gluon.
  pair<Energy,Energy> kplim = recoilPlusLimit(ph, pr, effMr);
  if ( kplim.first >= kplim.second ) return 0.0;

  // The lower limit is also the upper limit of the hard systems
  // positive momenum, which in turn, together with the gluons
  // momentum, gives an upper lower on the x-value for the extracted
  // parton. Veto if the PDF of the parton is zero.
  double xmin = 2.0*kplim.first/(rem->getBoost()*rem->parentMomentum()).plus();
  if ( rem->xfx(extracted, sqr(kt), xmin) <= 0.0 ) return 0.0;

  // The general suppression is alpha_S for the transverse momentum.
  double weight = Current<AriadneHandler>()->alphaS(sqr(kt));

  // This is the total incoming positive momentum.
  Energy pptot = ph.plus() + pr.plus();

  // This is the maximum positive momentum allowed for unsuppressed
  // emission from the incoming parton. If the whole available phase
  // space is unsuppressed no extra suppression is included.
  Energy kpsoftmax = pow(mu/kt, alpha)*pptot;
  if ( kpsoftmax >= kplim.second ) return 1.0;

  // We assume a flat distribution in positive momentum (ie. in
  // rapidity). Values larger than pptot(mu/kt)^alpha are suppressed
  // by a factor (ktmax/kt)^p where ktmax =
  // my(pptot/kp)^(1/alpha). The suppression is the ratio between the
  // integral of the suppressed and the unsuppressed phase space.
  if ( kpsoftmax <= kplim.first)
    weight *= (alpha/beta)*
      (pow(pptot/kplim.first, beta/alpha) -
       pow(pptot/kplim.second, beta/alpha))*pow(mu/kt, beta)/
      log(kplim.second/kplim.first);
  else
    weight *= ((alpha/beta)*(pow(pptot/kpsoftmax, beta/alpha) -
			     pow(pptot/kplim.second, beta/alpha))*
	       pow(mu/kt, beta) + log(kpsoftmax/kplim.first))/
      log(kplim.second/kplim.first);
  return weight;

}

pair<Energy,Energy> RemnantModel::
recoilPlusLimit(const LorentzMomentum & ph, const LorentzMomentum & pr,
		Energy effMr) const {
  // The limits are given such that the recoil gluon cannot have
  // apositive momentum larger than the remnant, or have a negative
  // momentum larger than the hard subsystem.
  pair<Energy,Energy> ret(Constants::MaxEnergy, ZERO);
  Energy2 kt2 = pr.perp2();
  Energy pptot = ph.plus() + pr.plus();
  Energy pmtot = ph.minus() + pr.minus();
  Energy2 mth2 = ph.mt2();
  Energy2 A = 2.0*kt2;
  Energy2 B = pptot*pmtot - mth2;
  Energy4 squarg = sqr(A + B) - 4.0*A*pptot*pmtot;
  if ( squarg <= ZERO ) return ret;
  Energy2 shmin = sqr(sqrt(pptot*pmtot) + effMr);
  Energy4 a = shmin*kt2;
  Energy2 b = shmin + kt2 - mth2;
  Energy4 squarg2 = sqr(b) - 4.0*a;
  if ( squarg2 <= ZERO ) return ret;
  
  ret.first = max((mth2 + kt2)/pmtot, 0.5*(b - sqrt(squarg2))/pmtot);
  ret.second = min((A + B + sqrt(squarg))/(4.0*pmtot),
		   0.5*(b + sqrt(squarg2))/pmtot);
  return ret;
}

Energy RemnantModel::effectiveMass(tcRemParPtr rem, tcPDPtr extracted) const {
  if ( !extracted ) extracted = &rem->extractedData();
  if ( rem->parentData().id() == extracted->id() ) return ZERO;
  return 2.0*rem->parentData().mass() + extracted->mass();
}

EmPtr RemnantModel::
generateRecoilGluon(const EmitterBase & em,
		    const QCDDipole & dip, tRemParPtr rem,
		    Energy rhomin, Energy rhomax) const {
 
  // First get the momentum of the remnant in the relevant system, and
  // check that it is inside limits.
  Lorentz5Momentum pr = rem->getBoost()*rem->momentum();
  Energy rho = pr.perp();
  if ( rho < rhomin ) return EmPtr();
  if ( rho > rhomax )
    Throw<InconsistentRecoilGluon>()
      << "The transverse momentum of one remnant was larger than the maximum "
      << "scale in the evolution. This is a serious errer. Consider contacting "
      << "the author." << Exception::runerror;

  RemnantGluonEmission e(em, dip);
  // There may be several dipoles generating this emission so we add a
  // little fluctuation in the scale to allow them to compete.
  e.rho = rho + rnd()*MeV;
  e.radiators.push_back(rem);
  e.colourParent = dip.iPart();
  e.antiColourParent = dip.oPart();
  e.mainParent = rem;
  e.pold = make_pair(rem->momentum(), rem->state()->hadronicMomentum());
  e.failsafe = true;
  return new_ptr(e);

}

bool RemnantModel::
performRecoilGluon(const RemnantGluonEmission & e) const {
  // First get the limits on the positive momentum of the recoil gluon
  // and check that they are consistent.
  tRemParPtr rem = dynamic_ptr_cast<tRemParPtr>(e.mainParent);
  LorentzRotation Rr = rem->getBoost();
  LorentzRotation Rir = rem->invBoost();
  Lorentz5Momentum pr = Rr*e.pold.first;
  LorentzMomentum ph = Rr*e.pold.second;
  Energy effMr = effectiveMass(rem);
  pair<Energy,Energy> pplim = recoilPlusLimit(ph, pr, effMr);
  double wmax =
    softSuppression(pr.perp(), pplim.first/(pr + ph).plus(), rem->x(),
		    rem->mu(), rem->alpha(), rem->beta());
  if ( pplim.first >= pplim.second || wmax <= 0.0 )
    Throw<InconsistentRecoilGluon>()
      << "The phase space for emitting a recoil gluon from a remnant with "
      << "non-zero transverse momentum was zero. This is a serious errer. "
      << "Consider contacting the author." << Exception::runerror;

  // Now generate the positive momentum of the recoil gluon.
  Energy plus = ZERO;
  LorentzMomentum pg;
  LorentzMomentum prnew;
  LorentzMomentum phnew;
  do {
    plus = pplim.first*pow(pplim.second/pplim.first, rnd());

    // Calculate the momentum of the gluon remnant and hard subsystem.
    pg = lightCone5(plus, pr.perp2()/plus, pr.x(), pr.y());
    Energy ppnew = (pr + ph).plus() - plus;
    Energy pmnew = (pr + ph).minus() - pg.minus();
    prnew = LorentzMomentum();
    phnew = ph;
    try {
      RemnantModel::setLightCones(prnew, phnew, ppnew, pmnew);
    }
    catch ( ImpossibleKinematics ) {
      cerr << "here we go" << endl;
      Throw<InconsistentRecoilGluon>()
	<< "The phase space for emitting a recoil gluon from a remnant with "
	<< "non-zero transverse momentum was zero. This is a serious errer. "
	<< "Consider contacting the author." << Exception::runerror;
    }

    double xnew = (plus + ph.plus())/(Rr*rem->parentMomentum()).plus();
    if ( rem->xfx(pr.perp2(), xnew) <= 0.0 ) continue;

  } while ( softSuppression(pr.perp(), plus/(pr + ph).plus(), rem->x(),
			    rem->mu(), rem->alpha(), rem->beta()) < rnd()*wmax );

  pg.transform(Rir);
  prnew.transform(Rir);
  phnew.transform(Rir);
  e.genmom = makeTriplet(prnew, pg, phnew);

  QCDDipole & d = dynamic_cast<QCDDipole &>(*e.dipole);
  e.Rh = rem->getHardTransform(e.pold.second, phnew);
  d.state()->transformHadronicState(e.Rh);
  tParPtr g = FSGluonEmitter::insertGluon(d, rem).first;

  g->momentum() = pg;
  rem->setMomentum(prnew);

  g->setVertex(rem->vertex());

  g->emission(&e);
  e.partons.push_back(g);

  return true;

}

void RemnantModel::
revertRecoilGluon(const RemnantGluonEmission & e) const {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*e.dipole);
  tRemParPtr rem = dynamic_ptr_cast<tRemParPtr>(e.mainParent);
  FSGluonEmitter::removeGluon(d, d.oPart(), rem);
  // LorentzMomentum phnew = e.genmom.third;
  // LorentzMomentum phold = e.pold.second;
  rem->setMomentum(e.pold.first);
  d.state()->transformHadronicState(e.Rh.inverse());
  d.state()->untouchHadronicState();
}

//
// Kinematics:
// 
// mt2 = xi*(1.0 - xi)*(mh2 + (1.0 - z)*Q2)/z + (mq2 - mh2)*(1.0 - xi)
//
// z = xi*(mh2 + Q2)*(1.0 - xi)/(mt2 + (mh2 - mq2)*(1.0 - xi) + Q2*xi*(1.0 -xi))
//
// -k2 = (mh2 + Q2)*(1.0 - xi)/z - mq2
//
// yq = log(mt2/(sqr(1.0 - xi)*s))/2.0
// yq = log(xi*((mh2 + (1.0 - z)*Q2)/z + (mq2 - mh2))/((1.0 - xi)*S))/2.0
//
// S = sqr(ph + pr)
//
// |d(mt2,yq)/d(xi,z)| =
//                mt2*(mh2 + Q2)/((1.0 - xi)*z*((mh2 + Q2)*(1.0 - z) + mq2*z))
//
// |d(mt2,yq)/d(k2,z)| =
//                mt2*(mh2 + Q2)/((k2 + mq2)*z*((mh2 + Q2)*(1.0 - z) + mq2*z))
//

EmPtr RemnantModel::
generateForcedHeavyQuark(const EmitterBase & em,
			 const QCDDipole & dip, tRemParPtr rem,
			 Energy rhomin, Energy rhomax) const {
  tcPDPtr q = rem->extractedData().CC();
  Energy mq = q->mass();
  Energy2 mq2 = sqr(mq);

  // Only force heavy quarks.
  if ( mq < Current<AriadneHandler>()->pTCut() ) return EmPtr();

  // If the remnant has a recoil pt larger than the quark mass, the
  // remnant should emit a recoil gluon first.
  if ( mq2 < rem->getPT2Kick() ) return EmPtr();

  if ( mq > rhomax )
    Throw<KinematicsException>()
      << "The forcing of an emission of a heavy quark from a remnant "
      "resulted in unordered scales. This should never happen. "
      "Consider contacting the author." << Exception::runerror;
  
  Lorentz5Momentum ph = rem->state()->hadronicMomentum();
  Lorentz5Momentum pr = rem->momentum();
  Energy2 mh2 = ph.mass2();
  double x = rem->x();
  tcPDPtr g = CurrentGenerator()->getParticleData(ParticleID::g);

  ISGtoQEmission e(em, dip, rem, q, mq, ph);
  Energy2 s = (pr + ph).m2();
  Energy2 mt2 = mq2;

  Energy2 Q2 = ZERO;
  if ( rem->state()->hadronicFS().size() == 1 ) {
    tRemParPtr hrem =
      dynamic_ptr_cast<tRemParPtr>(*(rem->state()->hadronicFS().begin()));
    if ( hrem && hrem->hard() ) Q2 = sqr(hrem->mu());
  }

  // Make sure that the scale is high enough for the PDF is non-zero.
  double maxpdfrat = 2.0*rem->xfratio(g, mt2, zmax(mt2, mq2, mh2, Q2));
  while ( maxpdfrat <= 0.0 ) {
    mt2 += sqr(Current<AriadneHandler>()->pTCut());
    maxpdfrat = 2.0*rem->xfratio(g, mt2, zmax(mt2, mq2, mh2, Q2));
  }

  while ( true ) {
    e.ymax = ymax(mq2, mq2, mh2, s);
    e.y = rnd(-e.ymax, e.ymax);
    e.xi = xiFn(mq2, e.y, s);
    e.z = zFn(mq2, e.xi, mq2, mh2, Q2);
    if ( mh2/e.xi + mq2/(1.0 - e.xi) >= mh2/x ) continue;
    double weight = (sqr(e.z) + sqr(1.0 - e.z));
    double pdfrat = min(rem->xfratio(g, mt2, e.z), maxpdfrat);
    if ( pdfrat <= 0.0 ) continue;
    weight *= pdfrat/maxpdfrat;
    // Jacobi determinant and reweighting.
    weight *= Jk2z(e.z, mq2, mh2, Q2);

    e.genmom = getMomenta(s, ZERO, e.xi, ZERO, mq2, mh2, pr, ph); 
 
    double recw = rem->recoilWeight(e.genmom.third, e.genmom.first, g);
    if ( recw <= 0.0 ) continue;
    recw /= rem->recoilWeight();
    if ( recw < 1.0 ) weight *= recw;

    if ( weight > 1.0 )
      Throw<WeightException>()
	<< "ISGtoQEmitter failed to overestimate the PDF ratio. "
	<< "If this hapens too often you should contact the author."
	<< Exception::warning;
      
    if ( weight < rnd() ) continue;

    e.ymax = log(sqrt(s)/Current<AriadneHandler>()->pTCut());
    e.rho = mq;
    e.radiators.push_back(rem);
    e.colourParent = e.mainParent = e.antiColourParent = rem;
    e.failsafe = true;
    return new_ptr(e);

  }

  return EmPtr();

}

bool RemnantModel::performGtoQ(const ISGtoQEmission & e) const {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*e.dipole);
  
  tRemParPtr rem = dynamic_ptr_cast<tRemParPtr>(e.mainParent);

  tcPDPtr g = CurrentGenerator()->getParticleData(ParticleID::g);
  setMom(d, e, rem, g);

  // Now create the new dipole.
  tQCDPtr newd = d.state()->create<QCDDipole>();
  newd->colourIndex(d.colourIndex());
  tParPtr q = d.state()->create(e.q, e.mainParent);
  q->momentum() = e.genmom.second;
  q->setVertex(e.mainParent->vertex());
  if ( rem == d.oPart() ) {
    newd->iPart(rem);
    newd->oPart(q);
    d.next(newd);
    newd->prev(&d);    
  } else {
    newd->oPart(rem);
    newd->iPart(q);
    d.prev(newd);
    newd->next(&d);    
  }
  newd->generateColourIndex();
  d.touch();
  rem->touch();
  q->emission(&e);
  e.partons.push_back(q);

  return true;

}

void RemnantModel::revertGtoQ(const ISGtoQEmission & e) const {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*e.dipole);
  tRemParPtr rem = dynamic_ptr_cast<tRemParPtr>(e.mainParent);

  if ( rem == d.oPart() ) {
    d.state()->forgetParton(d.next()->oPart());
    d.state()->forgetDipole(d.next());
    d.next(tQCDPtr());
  } else {
    d.state()->forgetParton(d.prev()->iPart());
    d.state()->forgetDipole(d.prev());
    d.prev(tQCDPtr());
  }
  d.untouch();
  rem->untouch();

  revertMom(d, e, rem);  
  d.state()->untouchHadronicState();
}

Triplet<Lorentz5Momentum,Lorentz5Momentum,Lorentz5Momentum>
RemnantModel::getMomenta(Energy2 S, Energy2 kt2, double xi, double phi,
			 Energy2 mq2, Energy2 mh2, const Lorentz5Momentum & pr,
			 const Lorentz5Momentum & ph) {
  Triplet<Lorentz5Momentum,Lorentz5Momentum,Lorentz5Momentum> p;
  Energy plus = sqrt(S);
  Energy pt = sqrt(kt2);
  Energy px = pt*cos(phi);
  Energy py = pt*sin(phi);
  p.third = Lorentz5Momentum(lightCone5((kt2 + mh2)/(plus*xi),
				     plus*xi, px, py), sqrt(mh2));
  p.second = Lorentz5Momentum(lightCone5((kt2 + mq2)/(plus*(1.0 - xi)),
				      plus*(1.0 - xi), -px, -py), sqrt(mq2));
  p.first = Lorentz5Momentum(lightCone5(plus - p.second.plus() - p.third.plus(),
				     0.0*GeV), ZERO);
  LorentzRotation R = Utilities::getBoostFromCM(make_pair(pr, ph));
  p.first.transform(R);
  p.second.transform(R);
  p.third.transform(R);
  
  return p;	    

}

void RemnantModel::
revertMom(QCDDipole & d, const ISQEmission & e, tRemParPtr rem) {
  // First set the new extracted parton and other kinematics.
  rem->setMomentum(e.pold.first, e.exorig);
  d.state()->transformHadronicState(e.Rh.inverse());
  rem->x(e.x);
}

void RemnantModel::
setMom(QCDDipole & d, const ISQEmission & e, tRemParPtr rem, tcPDPtr qex) {
  // First set the new extracted parton and other kinematics.
  e.Rh = rem->getHardTransform(e.pold.second, e.genmom.third);
  d.state()->transformHadronicState(e.Rh);
  rem->setMomentum(e.genmom.first, qex);
}




// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void RemnantModel::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void RemnantModel::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// The following static variable is needed for the type description
// system in ThePEG.
DescribeClass<RemnantModel,HandlerBase>
  describeRemnantModel("Ariadne5::RemnantModel", "libAriadne5.so");

void RemnantModel::Init() {

  static ClassDocumentation<RemnantModel> documentation
    ("RemnantModel is a helper class for emission models dealing with "
     "RemnantParton objects. This base class implements the default model "
     "for emitting a recoil gluon from a remnant that has received a "
     "transverse momentum kick. Furhter more it can emit initial-state "
     "type radiation between the hard sub-system and a remnant as well as "
     "vetoing or reweighting any other emission from a dipole connected to "
     "the remnant.");

}

