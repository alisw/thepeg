// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ISQtoGEmitter class.
//

#include "ISQtoGEmitter.h"
#include "ISQtoGEmission.h"
#include "RemnantModel.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "Ariadne/Cascade/AriadneHandler.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Current.h"
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

ISQtoGEmitter::ISQtoGEmitter() {}

ISQtoGEmitter::~ISQtoGEmitter() {}

IBPtr ISQtoGEmitter::clone() const {
  return new_ptr(*this);
}

IBPtr ISQtoGEmitter::fullclone() const {
  return new_ptr(*this);
}


bool ISQtoGEmitter::canHandle(const DipoleBase & dipole) const {
  tcQCDPtr d = dynamic_cast<const QCDDipole *>(&dipole);
  if ( !d ) return false;
  if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d->iPart()) )
    if ( !r->hard() && r->isG() ) return true;
  if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d->oPart()) )
    if ( !r->hard() && r->isG() ) return true;
  return false;
}

EmPtr ISQtoGEmitter::
generate(const DipoleBase & dipole, Energy rhomin, Energy rhomax) const {
  ISQtoGEmPtr esel;
  const QCDDipole & d = dynamic_cast<const QCDDipole &>(dipole);
  if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d.iPart()) )
    if ( !r->hard() && r->isG() ) esel = generate(r, d, rhomin, rhomax);
  if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d.oPart()) )
    if ( !r->hard() && r->isG() ) {
      if ( esel ) rhomin = max(rhomin, esel->rho);
      ISQtoGEmPtr e = generate(r, d, rhomin, rhomax);
      if ( !esel || ( e && e->rho > esel->rho ) ) esel = e;
    }

  return esel;

}

ISQtoGEmitter::ISQtoGEmPtr
ISQtoGEmitter::generate(tRemParPtr rem, const QCDDipole & d,
			Energy rhomin, Energy rhomax) const {
  static DebugItem debugweight("Ariadne5::ISQtoGEmitter::Weight");
  static DebugItem debugjac("Ariadne5::ISQtoGEmitter::Jacobian");
  ISQtoGEmPtr esel;
  // First collect some relevant data.
  Lorentz5Momentum ph = rem->state()->hadronicMomentum();
  Lorentz5Momentum pr = rem->momentum();
  Energy2 mh2 = ph.mass2();
  double x = rem->x();
  tcPDPtr g = &rem->extractedData();
  Energy2 s = (pr + ph).m2();

  ISQtoGEmission e(*this, d, rem, ph);

  const int nfl = Current<AriadneHandler>()->nFlav();

  for ( int ifl = 1; ifl <= nfl; ++ifl ) {

    e.q = generator()->getParticleData(rem == d.oPart()? -ifl: ifl);
    e.od = ( rem == d.iPart() )? d.prev(): d.next();
    e.mq = e.q->generateMass();
    Energy2 mq2 = sqr(e.mq);
    Energy2 mt2max = min(sqr(rhomax), sqr(s - mh2)/(4.0*s) + mq2);
    Energy2 mt2min = max(sqr(rhomin), sqr(handler().pTCut()) + mq2);
    mt2min = max(mt2min, sqr(e.rho));
    if ( mt2max <= mt2min ) continue;
    double C = 4.0/(3.0*Constants::pi);

    // We need to define the maximum possible value of the pdf
    // ratio. For efficiency, this may be done in different mt2 regions.
    PDFLimits limits = maxPDFRatios(rem, mt2max, mt2min, e.q, s, mh2, mq2);

    // Now we loop over regions.
    while ( limits.size() ) {
      mt2min = limits.back().first;
      double maxpdf = limits.back().second;
      limits.pop_back();
      Energy2 mt2 = mt2max;
      e.rho = ZERO;
      if ( maxpdf <= 0.0 ) continue;

      double yint = 2.0*RemnantModel::ymax(mt2min, mq2, mh2, s);
      double Jmargin = mt2min/(mt2min - mq2);

      while ( mt2 > mt2min ) {

	mt2 = rndsud(C*maxpdf*yint*Jmargin, mt2, mt2min);
	if ( mt2 < mt2min ) break;

	e.ymax = RemnantModel::ymax(mt2, mq2, mh2, s);
	e.y = 2.0*e.ymax*rnd() - e.ymax;
	e.xi = 1.0 - sqrt(mt2/s)*exp(-e.y);

	// Check if xi is kinematically possible
	if ( (mt2 + mh2 - mq2)/e.xi + mt2/(1.0 - e.xi) >= mh2/x ) continue;

	e.z = mh2*e.xi*(1.0 - e.xi)/(mh2*(1.0 - e.xi) + mq2*(e.xi - 1.0) + mt2);

	// Nominator of splitting function and correct PDF ratio.
	double pdfrat = rem->xfratio(e.q, mt2, e.z);
	if ( pdfrat <= 0.0 ) continue;
	if ( x > maxPDFX ) pdfrat = min(pdfrat, maxpdf);
	double weight = 0.5*(1.0 + sqr(1.0 - e.z))*pdfrat/maxpdf;

	// Jacobi determinant, correct for massless propagator and reweighting.
	double J = RemnantModel::Jk2z(e.z, mq2, mh2)/e.z;
	Energy2 v = (mh2*(1.0 - e.xi) + mt2 - mq2)/e.xi;
	J *= ((v + mq2)/v)/Jmargin;
	weight *= reweight(e)*J*2.0*e.ymax/yint;

	if ( debugweight && pdfrat > maxpdf ) {
	  cerr << "ISQtoGEmitter failed to overestimate PDF at " << endl
	       << " x = " << x << ", z = " << e.z << ", q = " << e.q->id()
	       << ", zmax = " << RemnantModel::zmax(mt2, mq2, mh2)
	       << ", mt = " << sqrt(mt2)/GeV
	       << ", ratio = " << pdfrat/maxpdf << endl;
	}
	if ( debugjac && J > 1.0 ) {
	  cerr << "ISQtoGEmitter found jacobian > 1 at " << endl
	       << " x = " << x << ", z = " << e.z << ", q = " << e.q->id()
	       << ", zmax = " << RemnantModel::zmax(mt2, mq2, mh2)
	       << ", mt = " << sqrt(mt2)/GeV << ", J = " << J << endl;
	}

	e.genmom = Current<AriadneHandler>()->remnantModel().getMomenta
	  (s, mt2 - mq2, e.xi, 2.0*Constants::pi*rnd(), mq2, mh2, pr, ph);

	if ( (rem->getBoost()*e.genmom.first).perp2() >= mt2 ) continue;
	double recw = rem->recoilWeight(e.genmom.third + e.genmom.second,
					e.genmom.first, e.q);
	if ( recw <= 0.0 ) continue;
	recw /= rem->recoilWeight();
	if ( recw < 1.0 ) weight *= recw;

	if ( weight > 1.0 )
	  Throw<WeightException>()
	    << "ISQtoGEmitter failed to overestimate the PDF ratio. "
	    << "If this hapens too often you should contact the author."
	    << Exception::warning;
      
	if ( weight < rnd() ) continue;

	e.rho = sqrt(mt2);
	e.ymax = log(sqrt(s)/handler().pTCut());
	e.radiators.push_back(rem);
	e.colourParent = e.mainParent = e.antiColourParent = rem;
	esel = new_ptr(e);
	break;

      }

      if ( e.rho > ZERO ) break;
      mt2max = mt2min;

    }
	
 
  }

  // *** ATTENTION *** Check this.

  return esel;

}

bool ISQtoGEmitter::
perform(const Emission & emission) const {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*emission.dipole);
  const ISQtoGEmission & e = dynamic_cast<const ISQtoGEmission &>(emission);

  tRemParPtr rem = dynamic_ptr_cast<tRemParPtr>(e.mainParent);

  if ( rem->recoilWeight(e.genmom.second + e.genmom.third,
			 e.genmom.first, e.q) <= 0.0 )
    Throw<WeightException>()
      << "Why wasn't this noticed in generate?." << Exception::runerror;
 
  setMom(d, e, rem, e.q);

  tParPtr q = d.state()->create(e.q, e.mainParent);
  q->momentum() = e.genmom.second;
  q->setVertex(e.mainParent->vertex());

  d.state()->sumTotalMomentum();
  if ( rem->recoilWeight() <= 0.0 )
    Throw<WeightException>()
      << "The transverse momentum of one remnant was larger than the maximum "
      << "scale in the evolution in ISQtoGEmitter. This is a serious errer. "
      << "Consider contacting the author." << Exception::runerror;
 

  if ( rem == d.iPart() ) {
    d.prev()->touch();
    d.prev()->next(tQCDPtr());
    d.prev(tQCDPtr());
    d.iPart(q);
  } else {
    d.next()->touch();
    d.next()->prev(tQCDPtr());
    d.next(tQCDPtr());
    d.oPart(q);
  }
  d.touch();
  rem->touch();
  q->emission(&e);
  e.partons.push_back(q);

  return true;

}

void ISQtoGEmitter::revert(const Emission & emission) const {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*emission.dipole);
  const ISQtoGEmission & e = dynamic_cast<const ISQtoGEmission &>(emission);
  tRemParPtr rem = dynamic_ptr_cast<tRemParPtr>(e.mainParent);

  if ( e.partons[0] == d.iPart() ) {
    d.prev(e.od);
    e.od->next(&d);
    d.iPart(rem);
  } else {
    d.next(e.od);
    e.od->prev(&d);
    d.oPart(rem);
  }
  d.state()->forgetParton(e.partons[0]);
  d.untouch();
  e.od->untouch();
  rem->untouch();

  revertMom(d, e, rem);
  d.state()->untouchHadronicState();

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ISQtoGEmitter::persistentOutput(PersistentOStream &) const {}

void ISQtoGEmitter::persistentInput(PersistentIStream &, int) {}

DescribeClass<ISQtoGEmitter,ISQEmitter>
describeAriadne5ISQtoGEmitter("Ariadne5::ISQtoGEmitter", "libAriadne5.so");

void ISQtoGEmitter::Init() {

  static ClassDocumentation<ISQtoGEmitter> documentation
    ("The ISQtoGEmitter class implements the initial-state splitting "
     "of a quark into an initial-state gluon and a final-state quark.");

}

