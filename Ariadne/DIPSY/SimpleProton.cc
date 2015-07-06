// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleProton class.
//

#include "SimpleProton.h"
#include "SimpleProtonState.h"
#include "DipoleEventHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/MaxCmp.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Utilities/Current.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "gsl/gsl_sf_erf.h"

using namespace DIPSY;

SimpleProton::SimpleProton()
  : theR(0.0*InvGeV), theR0(0.0*InvGeV), theRapidityWidth(0.0), theAngleWidth(0.0),
    theConnected(0), theCollapseTolerance(0.1*GeV) {}

SimpleProton::~SimpleProton() {}

IBPtr SimpleProton::clone() const {
  return new_ptr(*this);
}

IBPtr SimpleProton::fullclone() const {
  return new_ptr(*this);
}
InvEnergy SimpleProton::R() const {
  return theR > ZERO? theR: Current<DipoleEventHandler>()->rMax();
}

InvEnergy SimpleProton::r0() const {
  return theR0 > ZERO? theR0: Current<DipoleEventHandler>()->baryonSize();
}


void SimpleProton::initialize(const DipoleEventHandler & eh) {
  tcPDPtr gluon = generator()->getParticleData(ParticleID::g);
  if ( remdec && !remdec->canHandle(particle(), gluon ) )
    Throw<InitException>()
      << "The given remnant decayer could not handle the particle "
      << "in this wave function.";
}

Energy2 SimpleProton::m2() const {
  return sqr(particle()->mass());
}

DipoleStatePtr SimpleProton::
generate(const DipoleEventHandler & eh, Energy plus) {
  InvEnergy r = rndShiftGauss();
  Energy minus = m2()/plus;
  return new_ptr(SimpleProtonState(eh, plus, minus, theAngleWidth, theRapidityWidth,
				   connectOption(), new_ptr(WFInfo(this, r)), 1.0));
}

InvEnergy SimpleProton::rndShiftGauss() const {
  InvEnergy2 s1 = R()*r0()*sqrt(Constants::pi);
  InvEnergy2 s2 = sqr(R());
  InvEnergy2 s3 = s1*gsl_sf_erf(r0()/R());
  do {
    InvEnergy2 rsum = UseRandom::rnd()*(s1 + s2 + s3);
    if ( s1 > rsum )
      return R()*UseRandom::rndGauss() + r0();
    else if ( s1 + s2 > rsum )
      return R()*sqrt(-log(UseRandom::rnd())) + r0();
    InvEnergy r = r0() - R()*UseRandom::rndGauss();
    if ( r > 0.0*InvGeV && r/r0() > UseRandom::rnd() ) return r;
  } while ( true );
  return 0.0*InvGeV;
}

void SimpleProton::setParticle(PDPtr p) {
  if ( !BaryonMatcher::Check(abs(p->id())) )
    throw InterfaceException()
      << "Cannot set " << p->name()
      << " as particle for a proton wave function. "
      << "Only baryons are allowed." << Exception::warning;
  WaveFunction::setParticle(p);
}

void SimpleProton::fixValence(Step & step, tPPtr particle, const vector<PPtr> & val) const {

  if ( connectOption() != 0 || !remdec ) return;

  // If the hadron is not colour connected to the rest of the state,
  // then try to collapse it into the original hadron.
  if ( collapseToProton(step, particle, val) ) return;

  // Otherwise first find the gluon to split into a q and qq pair and split it.
  splitValence(step, particle, val);

}

bool SimpleProton::collapseToProton(Step & step, tPPtr p, const vector<PPtr> & val) const {
  if ( val.size() != 3 ) return false;
  if ( val[0]->colourLine() == val[1]->antiColourLine() ) {
    if ( val[1]->colourLine() != val[2]->antiColourLine() ) return false;
    if ( val[2]->colourLine() != val[0]->antiColourLine() ) return false;
  } else {
    if ( val[0]->colourLine() != val[2]->antiColourLine() ) return false;
    if ( val[2]->colourLine() != val[1]->antiColourLine() ) return false;
    if ( val[1]->colourLine() != val[0]->antiColourLine() ) return false;
  }

  // Get an instance of the hadron and check that its mass is not too different.
  Lorentz5Momentum pv = Utilities::sumMomentum(val);
  if ( abs(pv.mass() - p->mass()) > collapseTolerance() ) return false;

  // get the total momentum of the collision and determine how to
  // boost it to put the hadron on-shell.
  Lorentz5Momentum ptot = step.collision()->incoming().first->momentum() + 
                          step.collision()->incoming().second->momentum();
  Lorentz5Momentum pr = ptot - pv; // This is the recoil system to be used for momentum shuffling

  LorentzRotation Rshift;
  try {
    p->setMomentum(DipoleState::changeMass(pv, p->mass(), pr, &Rshift));
  } catch ( ImpossibleKinematics e ) {
    return false;
  }
  p->setVertex((val[0]->vertex() + val[1]->vertex() + val[2]->vertex())/3.0);
  step.addDecayProduct(val.begin(), val.end(), p);
  step.collision()->incoming().first->transform(Rshift);
  step.collision()->incoming().second->transform(Rshift);

  return true;

}

bool SimpleProton::splitValence(Step & step, tPPtr p, const vector<PPtr> & val) const {
  // First get the total momentum of the collision.
  LorentzMomentum ptot = step.collision()->incoming().first->momentum() + 
                          step.collision()->incoming().second->momentum();
  // Now we choose the valence gluon with the smallest transverse
  // momentum. Primarily we only consider those valence gluons which
  // can be split into massive quarks and di-quarks, if this is not
  // possible also others are considered with quark masses put to zero
  // and hoping for the best.
  MinCmp<Energy2,tPPtr> sel, primsel;
  LorentzPoint vnew;
  for ( int i = 0, N = val.size(); i < N; ++i ) {
    vnew += val[i]->vertex()/double(N);
    if ( val[i]->id() == ParticleID::g ) {
      sel(val[i]->momentum().perp2(), val[i]);
      if ( ptot.m() - (ptot -  val[i]->momentum()).m() >  2.0*p->mass() )
	primsel(val[i]->momentum().perp2(), val[i]);
    }
  }
  tPVector valence(1, sel.index());
  if ( valence.empty() ) return false;
  if ( primsel.index() ) valence[0] = primsel.index();

  tPPtr coln = valence[0]->colourNeighbour(val.begin(), val.end());
  tColinePtr coll = valence[0]->colourLine();
  tColinePtr acol = valence[0]->antiColourLine();
  if ( coln && (valence[0]->momentum() + coln->momentum()).m() < p->mass() ) {
    valence.push_back(coln);
    acol = coln->antiColourLine();
  }
  coln = valence[0]->antiColourNeighbour(val.begin(), val.end());
  if ( coln && (valence[0]->momentum() + coln->momentum()).m() < p->mass() ) {
    valence.push_back(coln);
    coll = coln->colourLine();
  }
  if ( valence.size() == 3 && ( valence[1] == valence[2] || coll == acol ) ) {
    valence = tPVector(1, valence[0]);
    coll = valence[0]->colourLine();
    acol = valence[0]->antiColourLine();
  }

  // Now find the momentum of the recoil system used to put the remnants on-shell.
  Lorentz5Momentum pv = Utilities::sumMomentum(valence);
  LorentzMomentum pr = ptot - pv;
  Energy maxmass = sqrt(sqr(ptot.m() - pr.mt()) - pr.perp2());

  const SimpleBaryonRemnantDecayer::BaryonContent & bi = remdec->getBaryonInfo(p->dataPtr());

  while ( true ) {
    // Select remnant flavours
    pair<int,int> r = bi.flavsel.select(UseRandom::current());
    PPtr qv = getParticleData(r.first*bi.sign)->produceParticle();
    PPtr dqv = getParticleData(r.second*bi.sign)->produceParticle();

    // Generate transverse momenum and energy splitting.
    TransverseMomentum ptv = remdec->pTGenerator()->generate();
    Energy2 mt02 = qv->momentum().mass2() + ptv.pt2();
    Energy2 mt12 = dqv->momentum().mass2() + ptv.pt2();
    // If still too large recoil mass we put quark masses to zero and hope for the best.
    if ( maxmass < 2.0*p->mass() ) mt02 = mt12 = ptv.pt2();
    double z = remdec->zGenerator().generate(dqv->dataPtr(), qv->dataPtr(), mt12);
    Energy mv = sqrt(mt02/(1.0 - z) + mt12/z);
    // Make shure the mass is small enough to be kinamatically allowd.
    if ( mv >= maxmass ) continue;

    // Calculate the momenta of the remnants and how much the rest
    // frame of the collision has moved due to momentum reshuffling.
    LorentzMomentum pqv = lightCone((1.0 - z)*mv, mt02/((1.0 - z)*mv), ptv);
    LorentzMomentum pdqv = lightCone(z*mv, mt12/(z*mv), -ptv);
    LorentzRotation Rshift;
    LorentzRotation Rr = Utilities::transformFromCMS(DipoleState::changeMass(pv, mv, pr, &Rshift));
    qv->setMomentum(Rr*pqv);
    if ( qv->hasColour() ) coll->addColoured(qv);
    else acol->addAntiColoured(qv);
    dqv->setMomentum(Rr*pdqv);
    if ( dqv->hasColour() ) coll->addColoured(dqv);
    else acol->addAntiColoured(dqv);
    qv->setVertex(vnew);
    dqv->setVertex(vnew);
    step.addDecayProduct(valence.begin(), valence.end(), qv);
    step.addDecayProduct(valence.begin(), valence.end(), dqv);
    step.collision()->incoming().first->transform(Rshift);
    step.collision()->incoming().second->transform(Rshift);
    return true;
  }

  return false;
}

void SimpleProton::persistentOutput(PersistentOStream & os) const {
  os << ounit(theR, InvGeV) << ounit(theR0, InvGeV)
     << ounit(theRapidityWidth, 1.0) << ounit(theAngleWidth, 1.0) << theConnected
     << ounit(theCollapseTolerance, GeV) << remdec;
}

void SimpleProton::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theR, InvGeV) >> iunit(theR0, InvGeV)
     >> iunit(theRapidityWidth, 1.0) >> iunit(theAngleWidth, 1.0) >> theConnected
     >> iunit(theCollapseTolerance, GeV) >> remdec;
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<SimpleProton,DIPSY::WaveFunction>
  describeDIPSYSimpleProton("DIPSY::SimpleProton", "libAriadne5.so libDIPSY.so");

void SimpleProton::Init() {

  static ClassDocumentation<SimpleProton> documentation
    ("The SimpleProton class represents the unevolved proton wave function "
     "described in terms of an equilatteral triangle of dipoles with a size "
     "distributed as a Gaussian.");

  static Parameter<SimpleProton,InvEnergy> interfaceR
    ("R",
     "The width of the Gaussian distribution in units of inverse GeV. If zero, "
     "the value of <interface>DipoleEventHandler::RMax</interface> of the "
     "controlling event handler will be used.",
     &SimpleProton::theR, InvGeV, 0.0*InvGeV, 0.0*InvGeV, 0.0*InvGeV,
     true, false, Interface::lowerlim);

  static Parameter<SimpleProton,InvEnergy> interfaceR0
    ("R0",
     "The shift in the average of the Gaussian distribution in units of "
     "inverse GeV. If zero, <interface>DipoleEventHandler::BaryonSize</interface> "
     "is used instead.",
     &SimpleProton::theR0, InvGeV, 0.0*InvGeV, 0.0*InvGeV, 0.0*InvGeV,
     true, false, Interface::lowerlim);

  static Parameter<SimpleProton, double> interfaceRapidityWidth
    ("RapidityWidth",
     "The width of the gaussian smearing in rapidity of the individual partons "
     "in the starting triangle.",
     &SimpleProton::theRapidityWidth, 1.0, 0.2, 0.0, 0.0,
     true, false, Interface::lowerlim);

  static Parameter<SimpleProton, double> interfaceAngleWidth
    ("AngleWidth",
     "The width in radians of the gaussian smearing of the shape of the triangle. ",
     &SimpleProton::theAngleWidth, 1.0, 0.2, 0.0, 0.0,
     true, false, Interface::lowerlim);

  static Switch<SimpleProton,int> interfaceConnected
    ("Connected",
     "If the three dipoles are connected or not.",
     &SimpleProton::theConnected, 0, true, false);
  static SwitchOption interfaceConnectedConnected
    (interfaceConnected,
     "Connected",
     "Three gluons.",
     0);
  static SwitchOption interfaceDisconnected
    (interfaceConnected,
     "Disconnected",
     "6 quarks, pairwise on top of each other.",
     1);

  static Parameter<SimpleProton,Energy> interfaceCollapseTolerance
    ("CollapseTolerance",
     "The maximum mass difference between the valence system and the "
     "hadron allowed for collapsing into the original hadron after "
     "the evolution.",
     &SimpleProton::theCollapseTolerance, GeV, 0.1*GeV, 0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

  static Reference<SimpleProton,SimpleBaryonRemnantDecayer> interfaceRemnantDecayer
    ("RemnantDecayer",
     "A RemnantDecayer object which is able to produce remnants for the particle ",
     &SimpleProton::remdec, false, false, true, true, true);

}

