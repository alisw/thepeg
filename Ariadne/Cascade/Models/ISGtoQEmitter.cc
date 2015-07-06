// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ISGtoQEmitter class.
//

#include "ISGtoQEmitter.h"
#include "ISGtoQEmission.h"
#include "RemnantModel.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "Ariadne/Cascade/AriadneHandler.h"
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

ISGtoQEmitter::ISGtoQEmitter() {}

ISGtoQEmitter::~ISGtoQEmitter() {}

IBPtr ISGtoQEmitter::clone() const {
  return new_ptr(*this);
}

IBPtr ISGtoQEmitter::fullclone() const {
  return new_ptr(*this);
}


bool ISGtoQEmitter::canHandle(const DipoleBase & dipole) const {
  tcQCDPtr d = dynamic_cast<const QCDDipole *>(&dipole);
  if ( !d ) return false;
  if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d->iPart()) )
    if ( !r->hard() && !r->isG() ) return true;
  if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d->oPart()) )
    if ( !r->hard() && !r->isG() ) return true;
  return false;
}

EmPtr ISGtoQEmitter::
generate(const DipoleBase & dipole, Energy rhomin, Energy rhomax) const {
  EmPtr esel;
  const QCDDipole & d = dynamic_cast<const QCDDipole &>(dipole);
  if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d.iPart()) )
    if ( !r->hard() && !r->isG() ) esel = generate(r, d, rhomin, rhomax);
  if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d.oPart()) )
    if ( !r->hard() && !r->isG() ) {
      if ( esel ) rhomin = max(rhomin, esel->rho);
      EmPtr e = generate(r, d, rhomin, rhomax);
      if ( !esel || ( e && e->rho > esel->rho ) ) esel = e;
    }

  return esel;

}

EmPtr ISGtoQEmitter::generate(tRemParPtr rem, const QCDDipole & d,
			      Energy rhomin, Energy rhomax) const {
  static DebugItem debugweight("Ariadne5::ISGtoQEmitter::Weight");
  static DebugItem debugjac("Ariadne5::ISGtoQEmitter::Jacobian");
  // First generated a possible forced emission if a heavy quark has
  // been emitted.
  const RemnantModel & remmod = Current<AriadneHandler>()->remnantModel();
  EmPtr esel = remmod.generateForcedHeavyQuark(*this, d, rem, rhomin, rhomax);
  if ( esel ) rhomin = max(rhomin, esel->rho);

  // First collect some relevant data.
  Lorentz5Momentum ph = rem->state()->hadronicMomentum();
  Lorentz5Momentum pr = rem->momentum();
  Energy2 mh2 = ph.mass2();
  double x = rem->x();
  tcPDPtr g = CurrentGenerator()->getParticleData(ParticleID::g);
  tcPDPtr q = rem->extractedData().CC();
  Energy mq = q->generateMass();
  Energy2 mq2 = sqr(mq);
  Energy2 s = (pr + ph).m2();
  Energy2 mt2max = min(sqr(rhomax), sqr(s - mh2)/(4.0*s) + mq2);
  if ( mt2max < mq2 )
    Throw<KinematicsException>()
      << "There was not enough energy for ISGtoQQEmitter to emit a quark. "
      << "This should never happen, execution will be aborted."
      << Exception::abortnow;
  Energy2 mt2min = max(sqr(rhomin), sqr(handler().pTCut()) + mq2);

  if ( mt2max <= mt2min ) return esel;

  ISGtoQEmission e(*this, d, rem, q, mq, ph);
  double C = 1.0/(4.0*Constants::pi);
  C *= preweight(e);


  // In the following we define z to be the ratio of the mass of the
  // hard subsystem before and after the emission. Also xi is defined
  // to be the fraction of the lightcone momenta of the incoming
  // particle taken by the original hard subsystem after the
  // emission. There are hard limits x<z<1 and x<xi<1.

  // We need to define the maximum possible value of the pdf
  // ratio. For efficiency, this may be done in different mt2 regions.
  PDFLimits limits = maxPDFRatios(rem, mt2max, mt2min, g, s, mh2, mq2);


  // Now we loop over regions.
  while ( limits.size() ) {

    mt2min = limits.back().first;
    // Increase overestimate at large x since this is a difficult region.
    double maxpdf = limits.back().second/(1.0 - x);
    limits.pop_back();
    Energy2 mt2 = mt2max;
    if ( maxpdf <= 0.0 ) continue;

    double yint = 2.0*RemnantModel::ymax(mt2min, mq2, mh2, s);

    while ( mt2 > mt2min ) {

      mt2 = rndsud(C*maxpdf*yint, mt2, mt2min);

      if ( mt2 < mt2min ) break;

      e.ymax = RemnantModel::ymax(mt2, mq2, mh2, s);
      e.y = 2.0*e.ymax*rnd() - e.ymax;

      e.xi = 1.0 - sqrt(mt2/s)*exp(-e.y);
      
      // Check if xi is kinematically possible
      if ( (mt2 + mh2 - mq2)/e.xi + mt2/(1.0 - e.xi) >= mh2/x ) continue;

      e.rho = sqrt(mt2);

      e.z = mh2*e.xi*(1.0 - e.xi)/(mh2*(1.0 - e.xi) + mq2*(e.xi - 1.0) + mt2);

      // Nominator of splitting function and correct PDF ratio.
      double pdfrat = rem->xfratio(g, mt2, e.z);
      if ( pdfrat <= 0.0 ) continue;
      if ( x > maxPDFX ) pdfrat = min(pdfrat, maxpdf);
      double weight = (sqr(e.z) + sqr(1.0 - e.z))*pdfrat/maxpdf;

      // Jacobi determinant and reweighting.
      double J = RemnantModel::Jk2z(e.z, mq2, mh2);
      weight *= reweight(e)*J*2.0*e.ymax/yint;

      if ( debugweight && pdfrat > maxpdf ) {
	cerr << "ISGtoQEmitter failed to overestimate PDF at " << endl
	     << " x = " << x << ", z = " << e.z << ", q = " << q->id()
 	     << ", zmax = " << RemnantModel::zmax(mt2, mq2, mh2)
	     << ", mt = " << sqrt(mt2)/GeV
	     << ", ratio = " << pdfrat/maxpdf << endl;
     }
      if ( debugjac && J > 1.0 ) {
	cerr << "ISGtoQEmitter found jacobian > 1 at " << endl
	     << " x = " << x << ", z = " << e.z << ", q = " << q->id()
 	     << ", zmax = " << RemnantModel::zmax(mt2, mq2, mh2)
	     << ", mt = " << sqrt(mt2)/GeV << ", J = " << J << endl;
      }

      e.genmom = remmod.getMomenta(s, mt2 - mq2, e.xi, 2.0*Constants::pi*rnd(),
				   mq2, mh2, pr, ph);

      if ( (rem->getBoost()*e.genmom.first).perp() >= e.rho ) continue;
      double recw =
	rem->recoilWeight(e.genmom.third + e.genmom.second, e.genmom.first, g);
      if ( recw <= 0.0 ) continue;
      recw /= rem->recoilWeight();
      if ( recw < 1.0 ) weight *= recw;

      if ( weight > 1.0 )
	Throw<WeightException>()
	  << "ISGtoQEmitter failed to overestimate the PDF ratio. "
	  << "If this hapens too often you should contact the author."
	  << Exception::warning;
      
      if ( weight < rnd() ) continue;

      e.ymax = log(sqrt(s)/handler().pTCut());
      e.radiators.push_back(rem);
      e.colourParent = e.mainParent = e.antiColourParent = rem;
      return new_ptr(e);

    }

    mt2max = mt2min;

  }

  // *** ATTENTION *** Check this.

  return esel;

}

bool ISGtoQEmitter::
perform(const Emission & emission) const {
  const ISGtoQEmission & e = dynamic_cast<const ISGtoQEmission &>(emission);
  return Current<AriadneHandler>()->remnantModel().performGtoQ(e);
}

void ISGtoQEmitter::revert(const Emission & emission) const {
  const ISGtoQEmission & e = dynamic_cast<const ISGtoQEmission &>(emission);
  Current<AriadneHandler>()->remnantModel().revertGtoQ(e);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ISGtoQEmitter::persistentOutput(PersistentOStream &) const {}

void ISGtoQEmitter::persistentInput(PersistentIStream &, int) {}

DescribeClass<ISGtoQEmitter,ISQEmitter>
describeAriadne5ISGtoQEmitter("Ariadne5::ISGtoQEmitter", "libAriadne5.so");

void ISGtoQEmitter::Init() {

  static ClassDocumentation<ISGtoQEmitter> documentation
    ("The ISGtoQEmitter class implements the initial-state splitting of gluons "
     "into a q-qbar pair.");

}

