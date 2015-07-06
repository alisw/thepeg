// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DYME class.
//

#include "DYME.h"
#include "RRGluonEmitter.h"
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

DYME::DYME() {}

DYME::~DYME() {}

IBPtr DYME::clone() const {
  return new_ptr(*this);
}

IBPtr DYME::fullclone() const {
  return new_ptr(*this);
}


bool DYME::canHandle(const DipoleBase & e) const {
  tcQCDPtr d = dynamic_cast<const QCDDipole *>(&e);
  if ( !d ) return false;
  tRemParPtr rem1 = dynamic_ptr_cast<tRemParPtr>(d->iPart());
  tRemParPtr rem3 = dynamic_ptr_cast<tRemParPtr>(d->oPart());
  if ( !rem1 || rem1->hard() || rem1->isG() || !rem1->untouched() ||
       !rem3 || rem3->hard() || rem3->isG() || !rem3->untouched() )
    return false;
  if ( d->state()->resonances().empty() ) return false;
  tcResPtr res = d->state()->resonances()[0];
  if ( res->orig()->id() != ParticleID::gamma &&
       res->orig()->id() != ParticleID::Z0 &&
       res->orig()->id() != ParticleID::Wplus &&
       res->orig()->id() != ParticleID::Wminus )
    return false;
  set<tPPtr> orig;
  for ( set<tParPtr>::const_iterator it = d->state()->hadronicFS().begin();
	it !=  d->state()->hadronicFS().end(); ++it )
    (**it).getOriginalParents(inserter(orig));
  if ( orig.empty() ) 
    Throw<Exception>()
      << "Ooops! DYME found no particles in the hard hadronic final state!"
      << Exception::abortnow;
  for ( set<tPPtr>::iterator it = orig.begin(); it != orig.end(); ++it )
    if ( !member(res->orig()->children(), *it) ) return false;

  return true;
}

bool DYME::overrides(const EmitterBase & em, DipoleBase &) const {
  if ( typeid(em) == typeid(FSGluonEmitter) ) return true;
  if ( typeid(em) == typeid(SpecialGluonEmitter) ) return true;
  if ( typeid(em) == typeid(RemnantGluonEmitter) ) return true;
  if ( typeid(em) == typeid(RRGluonEmitter) ) return true;
  if ( typeid(em) == typeid(ISGtoQEmitter) ) return true;
  return false;
}

EmPtr DYME::
generate(const DipoleBase & dipole, Energy rhomin, Energy rhomax) const {
  const QCDDipole & d = dynamic_cast<const QCDDipole &>(dipole);
  tRemParPtr rem1 = dynamic_ptr_cast<tRemParPtr>(d.iPart());
  tRemParPtr rem3 = dynamic_ptr_cast<tRemParPtr>(d.oPart());
  EmSel sel(rhomin);
  sel = ISGtoQEmitter::generate(rem1, d, rhomin, rhomax);
  sel = ISGtoQEmitter::generate(rem3, d, rhomin, rhomax);
  sel = generateG(rem1, rem3, d, rhomin, rhomax);

  return sel;

}

bool DYME::
perform(const Emission & emission) const {
  if ( tcISGtoQEmPtr e = dynamic_ptr_cast<tcISGtoQEmPtr>(&emission) )
    return Current<AriadneHandler>()->remnantModel().performGtoQ(*e);
  else
    return RRGluonEmitter::
      performRRG(dynamic_cast<const RRGluonEmission &>(emission));
}

void DYME::revert(const Emission & emission) const {
  if ( tcISGtoQEmPtr e = dynamic_ptr_cast<tcISGtoQEmPtr>(&emission) )
    Current<AriadneHandler>()->remnantModel().revertGtoQ(*e);
  else
    RRGluonEmitter::revertRRG(dynamic_cast<const RRGluonEmission &>(emission));
}


EmPtr DYME::generateG(tRemParPtr rem1, tRemParPtr rem3, const QCDDipole & d,
		       Energy rhomin, Energy rhomax) const {
  Energy2 S = d.sdip();
  if ( S <= sqr(2.0*rhomin) ) return EmPtr();
  Energy W = sqrt(S);
  RRGluonEmission e(*this, d, rem1, rem3, d.state()->hadronicMomentum());

  double C = preweight(e)*2.0/(3.0*Constants::pi);

  rhomax = min(rhomax, W/2.0);
  double yint = 2.0*acosh(0.5*W/rhomin);

  while (true) {
    if ( rhomax <= rhomin ) return EmPtr();
    double weight =
      RRGluonEmitter::gluonEngineRR(e, rhomin, rhomax, W, C, yint, 2, 2);
    if ( weight < 0.0 ) return EmPtr();
    rhomax = e.rho;
    if ( weight == 0.0 ) continue;
    weight *= reweight(e);

    // Do the actual matrix element correction
    Energy2 pt2 = sqr(e.rho);
    double yDY = rapidity(e.ophr);
    Energy2 shat = e.mh2 + 2.0*pt2 + 2.0*sqrt(pt2*(pt2 + e.mh2))*cosh(yDY - e.y);
    Energy2 that = -pt2 - sqrt(pt2*(pt2 + e.mh2))*exp(yDY - e.y);
    Energy2 uhat = e.mh2 - shat - that;
    double wme = sqr(e.rho)*(sqr(that) + sqr(uhat) + 2.0*shat*e.mh2)/
      (shat*that*uhat);
    weight *= wme/(sqr(e.x1) + sqr(e.x3));

    if ( weight > UseRandom::rnd() ) break;
  }

  e.colourParent = d.iPart();
  e.antiColourParent = d.oPart();
  e.ymax = log(W/Current<AriadneHandler>()->pTCut());
  e.radiators.push_back(rem1);
  e.radiators.push_back(rem3);

  return new_ptr(e);

}

double DYME::reweight(const Emission & emission) const {
  double w = 1.0;
  if ( const ISGtoQEmission * e =
       dynamic_cast<const ISGtoQEmission *>(&emission) ) {
    Energy2 mDY2 = e->pold.second.mass2();
    Energy2 mq2 = sqr(e->mq);
    Energy2 shat = mDY2/e->z;
    Energy2 that = -sqr(e->rho)/(1.0 - e->xi) + mq2;
    Energy2 uhat = mDY2 + mq2 - shat - that;
    w = 0.25*(sqr(shat) + sqr(uhat) + 2.0*mDY2*that)/
      (sqr(shat)*(sqr(e->z) + sqr(1.0 - e->z)));
    if ( w > 1.0 )
      cerr << "DYME gave matrix element correction above unity" << endl
	   << "mDY2 = " << mDY2/GeV2 << ", mq2 = " << mq2/GeV2
	   << ", shat = " << shat/GeV2 << ", that = " << that/GeV2
	   << ", uhat = " << uhat/GeV2
	   << ", mt2 = " << sqr(e->rho)/GeV2
	   << ", z = " << e->z
	   << ", xi = " << e->xi
	   << ", x = " << e->x << endl;
  }
  return w*ISGtoQEmitter::reweight(emission);
}

double DYME::preweight(const Emission & emission) const {
  double w = 1.0;
  if ( dynamic_cast<const ISGtoQEmission *>(&emission) ) w = 4.0;
  return w*ISGtoQEmitter::preweight(emission);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DYME::persistentOutput(PersistentOStream &) const {}

void DYME::persistentInput(PersistentIStream &, int) {}

DescribeClass<DYME,ISGtoQEmitter>
describeAriadne5DYME("Ariadne5::DYME", "libAriadne5.so");

void DYME::Init() {

  static ClassDocumentation<DYME> documentation
    ("The DYME class implements the emission of gluons and initial-state "
     "splitting of gluons into quarks according to the leading-order "
     "tree-level matrix element for deeply inelastic, neutral-current, "
     "lepton-proton scattering.");

}

