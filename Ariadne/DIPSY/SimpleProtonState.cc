// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleProtonState class.
//

#include "SimpleProtonState.h"
#include "ThePEG/Repository/UseRandom.h"
#include "SimpleProton.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Particle.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SimpleProtonState.tcc"
#endif


using namespace DIPSY;

SimpleProtonState::
SimpleProtonState(const DipoleEventHandler & eh, Energy plus, Energy minus,
		  double angleWidth, double rapWidth, int connected, WFInfoPtr wfi, double weight)
  : DipoleState(eh, wfi) {
  PPtr inc = wfi->wf().particle()->produceParticle(lightCone(plus, minus));
  inc->getInfo().push_back(wfi);
  vector<PartonPtr> & partons = theIncoming[inc];
  thePlus = plus;
  theMinus = minus;
  PartonPtr q1 = new_ptr(Parton());
  partons.push_back(q1);
  PartonPtr q2 = new_ptr(Parton());
  partons.push_back(q2);
  PartonPtr q3 = new_ptr(Parton());
  partons.push_back(q3);
  DipolePtr d1 = createDipole();
  DipolePtr d2 = createDipole();
  DipolePtr d3 = createDipole();
  d1->partons(make_pair(q1, q2));
  d2->partons(make_pair(q2, q3));
  d3->partons(make_pair(q3, q1));
  d1->neighbors(make_pair(d3, d2));
  d2->neighbors(make_pair(d1, d3));
  d3->neighbors(make_pair(d2, d1));
  generateColourIndex(d1);
  generateColourIndex(d2);
  generateColourIndex(d3);

  double phi1 = UseRandom::rnd()*2.0*Constants::pi;
  double phi2 = phi1 + 2.0*Constants::pi/3.0 + UseRandom::rndGauss(angleWidth,0.0);
  //  double phi3 = phi2 + 2.0*Constants::pi/3.0 + UseRandom::rndGauss(angleWidth,0.0);
  InvEnergy r = wfi->r()/sqrt(3.0);

  // gets p+ and y right. pT and p- adapts.
  //generate distribution in rapidity
  double R1(0.0),R2(0.0),R3(0.0);
  while (R1 > 1.0 || R1 == 0.0)
    R1 = exp(UseRandom::rndGauss(rapWidth,0.0))/3.0;
  while (R2 > (1.0 - R1) || R2 == 0.0)
    R2 = (1.0 - R1)*exp(UseRandom::rndGauss(rapWidth,0.0))/2.0;
  R3 = 1.0 - R1 - R2;
  q1->plus(plus*R1);
  q2->plus(plus*R2);
  q3->plus(plus*R3);

  q1->position(Parton::Point(r*cos(phi1), r*sin(phi1)));
  q2->position(Parton::Point(r*cos(phi2), r*sin(phi2)));
  //  q3->position(Parton::Point(r*cos(phi3), r*sin(phi3)));
  q3->position(-q1->position() - q2->position());

  q1->dipoles(make_pair(d3, d1));
  q2->dipoles(make_pair(d1, d2));
  q3->dipoles(make_pair(d2, d3));

  q1->pT(q1->recoil(q2) + q1->recoil(q3));
  q2->pT(q2->recoil(q1) + q2->recoil(q3));
  q3->pT(q3->recoil(q1) + q3->recoil(q2));

  double reScale = sqrt(minus/(sqr(q1->pT().pt())/q1->plus() +
			       sqr(q2->pT().pt())/q2->plus() + sqr(q3->pT().pt())/q3->plus()));

  q1->pT(q1->pT()*reScale);
  q2->pT(q2->pT()*reScale);
  q3->pT(q3->pT()*reScale);

  q1->y(log(q1->pT().pt()/q1->plus()));
  q2->y(log(q2->pT().pt()/q2->plus()));
  q3->y(log(q3->pT().pt()/q3->plus()));
  q1->oY(q1->y());
  q2->oY(q2->y());
  q3->oY(q3->y());
  q1->minus( q1->pT().pt()*exp( q1->y() ) );
  q2->minus( q2->pT().pt()*exp( q2->y() ) );
  q3->minus( q3->pT().pt()*exp( q3->y() ) );

  q1->valence(true);
  q2->valence(true);
  q3->valence(true);

  q1->valencePT(q1->pT());
  q2->valencePT(q2->pT());
  q3->valencePT(q3->pT());
  q1->valencePlus(q1->plus());
  q2->valencePlus(q2->plus());
  q3->valencePlus(q3->plus());
  theInitialDipoles.push_back(d1);
  theInitialDipoles.push_back(d2);
  theInitialDipoles.push_back(d3);

  if ( connected == 1 ) {//disconnected
    PartonPtr qbar1 = new_ptr(Parton());
    partons.push_back(qbar1);
    PartonPtr qbar2 = new_ptr(Parton());
    partons.push_back(qbar2);
    PartonPtr qbar3 = new_ptr(Parton());
    partons.push_back(qbar3);
    d1->partons(make_pair(q1, qbar2));
    d2->partons(make_pair(q2, qbar3));
    d3->partons(make_pair(q3, qbar1));
    d1->neighbors(make_pair(tDipolePtr(), tDipolePtr()));
    d2->neighbors(make_pair(tDipolePtr(), tDipolePtr()));
    d3->neighbors(make_pair(tDipolePtr(), tDipolePtr()));

    q1->plus(q1->plus()/2.0);
    q2->plus(q2->plus()/2.0);
    q3->plus(q3->plus()/2.0);
    qbar1->plus(q1->plus());
    qbar2->plus(q2->plus());
    qbar3->plus(q3->plus());

    Parton::Point pos1 = Parton::Point(q1->position().x() + UseRandom::rndGauss(0.0,0.0)/GeV,
				       q1->position().y() + UseRandom::rndGauss(0.0,0.0)/GeV);
    Parton::Point pos2 = Parton::Point(q2->position().x() + UseRandom::rndGauss(0.0,0.0)/GeV,
				       q2->position().y() + UseRandom::rndGauss(0.0,0.0)/GeV);
    Parton::Point pos3 = Parton::Point(q3->position().x() + UseRandom::rndGauss(0.0,0.0)/GeV,
				       q3->position().y() + UseRandom::rndGauss(0.0,0.0)/GeV);

    qbar1->position(pos1);
    qbar2->position(pos2);
    qbar3->position(pos3);
    qbar1->dipoles(make_pair(d3, tDipolePtr()));
    qbar2->dipoles(make_pair(d1, tDipolePtr()));
    qbar3->dipoles(make_pair(d2, tDipolePtr()));
    q1->dipoles(make_pair(tDipolePtr(), d1));
    q2->dipoles(make_pair(tDipolePtr(), d2));
    q3->dipoles(make_pair(tDipolePtr(), d3));
    q1->pT(q1->pT()/2.0);
    q2->pT(q2->pT()/2.0);
    q3->pT(q3->pT()/2.0);
    qbar1->pT(q1->pT()/2.0);
    qbar2->pT(q2->pT()/2.0);
    qbar3->pT(q3->pT()/2.0);
    q1->updateYMinus();
    q2->updateYMinus();
    q3->updateYMinus();
    qbar1->updateYMinus();
    qbar2->updateYMinus();
    qbar3->updateYMinus();
    qbar1->oY(qbar1->y());
    qbar2->oY(qbar2->y());
    qbar3->oY(qbar3->y());

    qbar1->valence(true);
    qbar2->valence(true);
    qbar3->valence(true);

    q1->valencePT(q1->pT());
    q2->valencePT(q2->pT());
    q3->valencePT(q3->pT());
    qbar1->valencePT(qbar1->pT());
    qbar2->valencePT(qbar2->pT());
    qbar3->valencePT(qbar3->pT());
    q1->valencePlus(q1->plus());
    q2->valencePlus(q2->plus());
    q3->valencePlus(q3->plus());
    qbar1->valencePlus(qbar1->plus());
    qbar2->valencePlus(qbar2->plus());
    qbar3->valencePlus(qbar3->plus());

    qbar1->flavour(ParticleID::dbar);
    q1->flavour(ParticleID::d);
    qbar2->flavour(ParticleID::ubar);
    q2->flavour(ParticleID::u);
    qbar3->flavour(ParticleID::ubar);
    q3->flavour(ParticleID::u);
    if ( wfi->wf().particle() ) {
      if ( wfi->wf().particle()->id() == ParticleID::pplus )
	q1->flavour(ParticleID::u);
      if ( wfi->wf().particle()->id() == ParticleID::pbarminus )
	qbar1->flavour(ParticleID::ubar);
    }
  }
}

SimpleProtonState::~SimpleProtonState() {}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeNoPIOClass<SimpleProtonState,DIPSY::DipoleState>
  describeDIPSYSimpleProtonState("DIPSY::SimpleProtonState", "libAriadne5.so libDIPSY.so");


void SimpleProtonState::Init() {}

