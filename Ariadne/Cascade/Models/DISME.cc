// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISME class.
//

#include "FSGluonEmission.h"
#include "DISME.h"
#include "FSGluonEmitter.h"
#include "SpecialGluonEmitter.h"
#include "RemnantGluonEmitter.h"
#include "RemnantModel.h"
#include "ISGtoQEmitter.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "Ariadne/Cascade/AriadneHandler.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/DebugItem.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

DISME::DISME() {}

DISME::~DISME() {}

IBPtr DISME::clone() const {
  return new_ptr(*this);
}

IBPtr DISME::fullclone() const {
  return new_ptr(*this);
}


bool DISME::canHandle(const DipoleBase & e) const {
  tcQCDPtr d = dynamic_cast<const QCDDipole *>(&e);
  if ( !d ) return false;
  if ( d->state()->hadronicFS().size() != 1 ) return false;
  tRemParPtr rhard = dynamic_ptr_cast<tRemParPtr>(d->iPart());
  tRemParPtr rsoft = dynamic_ptr_cast<tRemParPtr>(d->oPart());
  if ( !rhard || ! rsoft ) return false;
  if ( !rhard->hard() ) swap(rhard, rsoft);
  if ( !rhard->hard() ) return false;
  if ( rsoft->hard() || !rsoft->untouched() ) return false;
  if ( *(d->state()->hadronicFS().begin()) != rhard ) return false;
  return true;
}

bool DISME::overrides(const EmitterBase & em, DipoleBase &) const {
  if ( typeid(em) == typeid(FSGluonEmitter) ) return true;
  if ( typeid(em) == typeid(SpecialGluonEmitter) ) return true;
  if ( typeid(em) == typeid(RemnantGluonEmitter) ) return true;
  if ( typeid(em) == typeid(ISGtoQEmitter) ) return true;
  return false;
}

EmPtr DISME::
generate(const DipoleBase & dipole, Energy rhomin, Energy rhomax) const {
  const QCDDipole & d = dynamic_cast<const QCDDipole &>(dipole);
  tRemParPtr rhard = dynamic_ptr_cast<tRemParPtr>(d.iPart());
  tRemParPtr rsoft = dynamic_ptr_cast<tRemParPtr>(d.oPart());
  if ( !rhard->hard() ) swap(rhard, rsoft);
  EmSel sel(rhomin);
  sel = generateG(rsoft, rhard, d, rhomin, rhomax);
  sel = generateQ(rsoft, rhard, d, rhomin, rhomax);
  return sel;

  return EmPtr();

}


bool DISME::
perform(const Emission & emission) const {
  if ( tcRemGEmPtr e = dynamic_ptr_cast<tcRemGEmPtr>(&emission) )
    return performG(*e);
  if ( tcISGtoQEmPtr e = dynamic_ptr_cast<tcISGtoQEmPtr>(&emission) )
    return performQ(*e);

  return false;
}

void DISME::revert(const Emission & emission) const {
  if ( tcRemGEmPtr e = dynamic_ptr_cast<tcRemGEmPtr>(&emission) )
    revertG(*e);
  if ( tcISGtoQEmPtr e = dynamic_ptr_cast<tcISGtoQEmPtr>(&emission) )
    revertQ(*e);
}


EmPtr DISME::generateG(tRemParPtr rsoft, tRemParPtr rhard, const QCDDipole & d,
		       Energy rhomin, Energy rhomax) const {
  DebugItem tuple("Ariadne5::DISME::Tuple");
  LorentzRotation Rr = rsoft->getBoost();
  Lorentz5Momentum ph = rhard->momentum();
  Lorentz5Momentum pr = rsoft->momentum();
  Lorentz5Momentum phr = Rr*rhard->momentum();
  Lorentz5Momentum prr = Rr*rsoft->momentum();
  Energy2 S = (ph + pr).m2();
  Energy2 mq2 = ph.mass2();
  Energy W = sqrt(S);
  Energy2 Q2 = sqr(rhard->mu());
  double nu = rhard->y();
  RemnantGluonEmission e(*this, d, prr.mass2()/S, mq2/S);

  double C = preweight(e)*2.0/(3.0*Constants::pi);
  rhomax = min(rhomax, W/2.0);
  double yint = 2.0*acosh(0.5*W/rhomin);

  while (true) {
    double weight =
      FSGluonEmitter::gluonEngine(e, rhomin, rhomax, W, C, yint, 2, 2);
    if ( weight < 0.0 ) return EmPtr();
    rhomax = e.rho;
    if ( weight == 0.0 ) continue;
    weight *= reweight(e);

    // Check if this is to be treated as a final-state splitting.
    //    bool isFS = e.rho*exp(e.y) < Q2/W;
    bool isFS = true;

    try {
      
      e.genmom =
	getMomenta(S, e.x1, e.x3, prr.mass(), ZERO, sqrt(mq2), true, !isFS);
      if ( e.genmom.second.plus() > Q2/W ) isFS = false;

      e.genmom =
	getMomenta(S, e.x1, e.x3, prr.mass(), ZERO, sqrt(mq2),
		   true, !isFS, UseRandom::rnd(2.0*Constants::pi), pr, ph);
    
      // Mass of the hard system after the gluon is emitted.
      Energy ppq = (mq2 + sqr(e.rho))/(W - e.rho*exp(-e.y));
      Energy2 sh = W*(ppq + e.rho*exp(e.y));

      // Check if there is enough phase space available.
      if ( sh >= S ) continue;

      double z = (mq2 + Q2)/(sh + Q2);
      double xi = 1.0 - e.rho*exp(-e.y)/W;

      // Matrix element multiplied by (1-z)(1-xi)/2
      weight *= (sqr(z) + sqr(xi) +
		 (2.0 + xi*z*(2.0 + 8.0*(1.0 - nu)/(1.0 + sqr(1.0 - nu))))*
		 (1.0 - xi)*(1.0 - z))/2.0;
      // Jacobi determinant divided by z(1-z)(1-xi)
      weight *= S/(S + mq2*exp(-2.0*e.y));

      LorentzMomentum pem = Rr*e.genmom.second;
      double remw = rsoft->softSuppression(e.rho, pem.plus()/(prr + phr).plus());
      remw = min(remw,
		 rhard->softSuppression(e.rho, pem.minus()/(prr + phr).minus()));
      // If the remnant has acquired a transverse momentum an extra
      // suppression is needed for the subsequent emission of a recoil
      // gluon.  The returned weight is the ratio of the new and the old
      // suppression But also if the transverse momentum of the recoil
      // gluon is above the scale of the emission we must veto.
      if ( (Rr*e.genmom.first).perp() >= e.rho ) continue;
      double recw =
	rsoft->recoilWeight(e.genmom.second + e.genmom.third, e.genmom.first);
      if ( recw <= 0.0 ) continue;
      recw /= rsoft->recoilWeight();
      if ( recw < 1.0 ) remw *= recw;
      if ( !isFS ) weight *= remw;

    } catch ( ImpossibleKinematics ) {
      continue;
    }

    if ( weight < UseRandom::rnd() )  continue;

    e.mainParent = isFS? rhard: rsoft;
    e.colourParent = d.iPart();
    e.antiColourParent = d.oPart();
    e.pold = make_pair(rsoft->momentum(), rhard->momentum());
    e.ymax = log(W/Current<AriadneHandler>()->pTCut());
    e.radiators.push_back(rsoft);
    e.radiators.push_back(rhard);
    if ( tuple )
      cerr << e.rho/GeV << " " << e.y << " "
	   << sqrt(Q2)/GeV << " " << W/GeV << endl;

    return new_ptr(e);

  }

  return EmPtr();

}

bool
DISME::performG(const RemnantGluonEmission & e) const {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*e.dipole);
  tRemParPtr rsoft = dynamic_ptr_cast<tRemParPtr>(e.radiators[0]);
  tRemParPtr rhard = dynamic_ptr_cast<tRemParPtr>(e.radiators[1]);
  tParPtr g = FSGluonEmitter::insertGluon(d, e.mainParent).first;
  rsoft->setMomentum(e.genmom.first);
  g->momentum() = e.genmom.second;
  rhard->setMomentum(e.genmom.third);

  g->setVertex(e.mainParent->vertex());

  g->emission(&e);
  e.partons.push_back(g);

  return true;  

}

void DISME::revertG(const RemnantGluonEmission & e) const {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*e.dipole);
  tRemParPtr rsoft = dynamic_ptr_cast<tRemParPtr>(e.radiators[0]);
  tRemParPtr rhard = dynamic_ptr_cast<tRemParPtr>(e.radiators[1]);
  FSGluonEmitter::removeGluon(d, d.oPart(), e.mainParent);
  rsoft->setMomentum(e.pold.first);
  rhard->setMomentum(e.pold.second);
  d.state()->untouchHadronicState();
}


EmPtr DISME::generateQ(tRemParPtr rsoft, tRemParPtr rhard, const QCDDipole & d,
		       Energy rhomin, Energy rhomax) const {
  static DebugItem debugweight("Ariadne5::DISME::Weight");
  static DebugItem debugjac("Ariadne5::DISME::Jacobian");
  const RemnantModel & remmod = Current<AriadneHandler>()->remnantModel();
  // First generated a possible forced emission if a heavy quark has
  // been emitted.
  EmPtr esel = Current<AriadneHandler>()->remnantModel().
    generateForcedHeavyQuark(*this, d, rsoft, rhomin, rhomax);
  if ( esel ) rhomin = max(rhomin, esel->rho);

  Lorentz5Momentum ph = rhard->momentum();
  Lorentz5Momentum pr = rsoft->momentum();
  Energy2 S = (ph + pr).m2();
  Energy2 mh2 = ph.mass2();
  Energy2 Q2 = sqr(rhard->mu());
  double nu = rhard->y();
  tcPDPtr g = CurrentGenerator()->getParticleData(ParticleID::g);
  tcPDPtr q = &rsoft->data();
  Energy mq = q->generateMass();
  Energy2 mq2 = sqr(mq);
  double x = rsoft->x();

  Energy2 mt2max = min(sqr(rhomax), sqr(S + mh2 - mq2)/(4.0*S));
  if ( mt2max <= mq2 )
    Throw<KinematicsException>()
      << "There was not emough energy for ISGtoQQEmitter to emit a quark. "
      << "This should never happen, execution will be aborted."
      << Exception::abortnow;
  Energy2 mt2min = max(sqr(rhomin), sqr(handler().pTCut()) + mq2);

  if ( mt2max <= mt2min ) return esel;

  ISGtoQEmission e(*this, d, rsoft, q, mq, ph);

  double C = preweight(e)/(4.0*Constants::pi);
  // In the following we define z to be the ratio of the mass of the
  // hard subsystem before and after the emission. Also xi is defined
  // to be the fraction of the lightcone momenta of the incoming
  // particle taken by the original hard subsystem after the
  // emission. There are hard limits x<z<1 and x<xi<1.

  // We need to define the maximum possible value of the pdf
  // ratio. For efficiency, this may be done in different mt2 regions.
  PDFLimits limits = maxPDFRatios(rsoft, mt2max, mt2min, g, S, mh2, mq2, Q2);

  // Now we loop over regions.
  while ( limits.size() ) {

    mt2min = limits.back().first;
    double maxpdf = limits.back().second;
    limits.pop_back();
    Energy2 mt2 = mt2max;
    double yint = 2.0*RemnantModel::ymax(mt2min, mq2, mh2, S);
    if ( maxpdf <= 0.0 ) continue;

    while ( true ) {

      mt2 = rndsud(C*maxpdf, mt2, mt2min);

      if ( mt2 < mt2min ) break;

      e.ymax = RemnantModel::ymax(mt2, mq2, mh2, S);
      e.y = 2.0*e.ymax*rnd() - e.ymax;

      //      Energy2 sh = mt2*sqr(2*cosh(e.y));
      // check if y and pt2 is kinematically possible
      //      if ( sh > S ) continue;

      //      e.z = (mq2 + Q2)/(sh + Q2);
      //      e.xi = 1.0/(1.0 + exp(2.0*e.y));
      e.xi = 1.0 - sqrt(mt2/S)*exp(-e.y);
      e.z = (mh2 + Q2)*e.xi*(1.0 - e.xi)/
      	(mt2 + (mh2 - mq2)*(1.0 - e.xi) + Q2*e.xi*(1.0 - e.xi));

      Energy2 sh = ((1.0 - e.xi)*mt2 + mh2 - mq2)/e.xi;
      if ( sh >= S ) continue;

      double weight = (sqr(e.z) + sqr(1.0 - e.z))*e.xi/(1.0 - e.xi) +
	8.0*e.z*(1.0 - e.z)*(1.0 - nu)/(1.0 + sqr(1.0 + nu));
      weight *= reweight(e)*2.0*e.ymax/yint;
      double pdfrat = rsoft->xfratio(g, mt2, e.z);
      if ( pdfrat <= 0.0 ) continue;
      if ( x > maxPDFX ) pdfrat = min(pdfrat, maxpdf);
      weight *= pdfrat/maxpdf;

      // multiply by jacobi determinant / 2.0
      //      weight *= e.z*mt2/(sh + Q2);
      double J = RemnantModel::Jxiz(e.xi, e.z, mq2, mh2, Q2);
      weight *= J;

      if ( debugweight && pdfrat > maxpdf ) {
	cerr << "DISME failed to overestimate PDF at " << endl
	     << " x = " << x << ", z = " << e.z << ", q = " << q->id()
 	     << ", zmax = " << RemnantModel::zmax(mt2, mq2, mh2, Q2)
	     << ", mt = " << sqrt(mt2)/GeV
	     << ", ratio = " << pdfrat/maxpdf << endl;
      }
      if ( debugjac && J > 1.0 - e.xi ) {
	cerr << "DISME found jacobian > 1 at " << endl
	     << " x = " << x << ", z = " << e.z << ", q = " << q->id()
 	     << ", zmax = " << RemnantModel::zmax(mt2, mq2, mh2, Q2)
	     << ", mt = " << sqrt(mt2)/GeV << ", J = " << J << endl;
      }

      // We now place the emission in the dipole cms as if the remnant
      // didn't have any pt before, and will not acquire any after the
      // emission.
      e.genmom = remmod.getMomenta(S, mt2 - mq2, e.xi, 2.0*Constants::pi*rnd(),
				   mq2, mh2, pr, ph);

      e.rho = sqrt(mt2);
      if ( (rsoft->getBoost()*e.genmom.first).perp() >= e.rho ) continue;
      double recw =
	rsoft->recoilWeight(e.genmom.third + e.genmom.second, e.genmom.first, g);
      if ( recw <= 0.0 ) continue;
      recw /= rsoft->recoilWeight();
      if ( recw < 1.0 ) weight *= recw;



      if ( weight > 1.0 )
	Throw<WeightException>()
	  << "DISME failed to overestimate the PDF ratio. "
	  << "If this hapens too often you should contact the author."
	  << Exception::warning;
      
      if ( weight < rnd() ) continue;

      e.pold = make_pair(rsoft->momentum(), rhard->momentum());
      e.radiators.push_back(rsoft);
      e.radiators.push_back(rhard);
      e.mainParent = rsoft;
      e.colourParent = d.iPart();
      e.antiColourParent = d.oPart();
      return new_ptr(e);

    }

    mt2max = mt2min;

  }

  // *** ATTENTION *** Check this.
  return esel;

}

bool
DISME::performQ(const ISGtoQEmission & e) const {
  if ( e.radiators.size() < 2 )
    return Current<AriadneHandler>()->remnantModel().performGtoQ(e);
  QCDDipole & d = dynamic_cast<QCDDipole &>(*e.dipole);
  tRemParPtr rsoft = dynamic_ptr_cast<tRemParPtr>(e.mainParent);
  tRemParPtr rhard = dynamic_ptr_cast<tRemParPtr>(e.radiators[1]);

  tcPDPtr g = CurrentGenerator()->getParticleData(ParticleID::g);
  rsoft->setMomentum(e.genmom.first, g);
  rhard->setMomentum(e.genmom.third);

  // Now create the new dipole.
  tQCDPtr newd = d.state()->create<QCDDipole>();
  newd->colourIndex(d.colourIndex());
  tParPtr q = d.state()->create(e.q, e.mainParent);
  q->momentum() = e.genmom.second;
  q->setVertex(e.mainParent->vertex());
  if ( rsoft == d.oPart() ) {
    newd->iPart(rsoft);
    newd->oPart(q);
    d.next(newd);
    newd->prev(&d);    
  } else {
    newd->oPart(rsoft);
    newd->iPart(q);
    d.prev(newd);
    newd->next(&d);    
  }
  newd->generateColourIndex();
  d.touch();
  rsoft->touch();
  rhard->touch();
  q->emission(&e);
  e.partons.push_back(q);

  return true;

}

void DISME::revertQ(const ISGtoQEmission & e) const {
  if ( e.radiators.size() < 2 )
    return Current<AriadneHandler>()->remnantModel().revertGtoQ(e);
  QCDDipole & d = dynamic_cast<QCDDipole &>(*e.dipole);
  tRemParPtr rsoft = dynamic_ptr_cast<tRemParPtr>(e.mainParent);
  tRemParPtr rhard = dynamic_ptr_cast<tRemParPtr>(e.radiators[1]);

  if ( rsoft == d.oPart() ) {
    d.state()->forgetParton(d.next()->oPart());
    d.state()->forgetDipole(d.next());
    d.next(tQCDPtr());
  } else {
    d.state()->forgetParton(d.prev()->iPart());
    d.state()->forgetDipole(d.prev());
    d.prev(tQCDPtr());
  }
  rsoft->setMomentum(e.pold.first, e.exorig);
  rhard->setMomentum(e.pold.second);
  rsoft->x(e.x);

  d.untouch();
  rsoft->untouch();
  rhard->untouch();
  d.state()->untouchHadronicState();
  
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DISME::persistentOutput(PersistentOStream &) const {}

void DISME::persistentInput(PersistentIStream &, int) {}

DescribeClass<DISME,ISQEmitter>
describeAriadne5DISME("Ariadne5::DISME", "libAriadne5.so");

void DISME::Init() {

  static ClassDocumentation<DISME> documentation
    ("The DISME class implements the emission of gluons and initial-state "
     "splitting of gluons into quarks according to the leading-order "
     "tree-level matrix element for deeply inelastic, neutral-current, "
     "lepton-proton scattering.");

}

