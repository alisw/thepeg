// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TestDecays class.
//

#include "TestDecays.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TestDecays.tcc"
#endif


using namespace ThePEG;

TestDecays::~TestDecays() {}

NoPIOClassDescription<TestDecays> TestDecays::initTestDecays;
// Definition of the static class description member.

Energy2 ktclus2(const LorentzMomentum & pi, const LorentzMomentum & pj) {
  double deta2 = sqr(pi.eta() - pj.eta());
  double dphi = abs(pi.phi() - pj.phi());
  if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
  return min(pi.perp2(), pj.perp2())*(deta2 + sqr(dphi));
}

Energy2 durham(const LorentzMomentum & pi, const LorentzMomentum & pj) {
  return 2.0*sqr(min(pi.e(), pj.e()))*(1.0 - cos(pi.angle(pj)));
}

void TestDecays::Init() {
  breakThePEG();
  if ( N() <= 0 ) N(10000);
  //  testGG(eventGenerator(), 3*GeV, N());
  testDecay(eventGenerator(), ParticleID::D_splus, N());
//   testDecay(eventGenerator(), ParticleID::Bbar0, N());
//   testDecay(eventGenerator(), ParticleID::Lambda_cplus, N());
//   testDecay(eventGenerator(), ParticleID::Lambda_cbarminus, N());
//   testDecay(eventGenerator(), ParticleID::Upsilon, N());
/*  ofstream massive("massive.dat");
  ofstream massless("massless.dat");
  UseRandom rnd(&(eventGenerator()->random()));
  int np = 7;
  for ( int iev = 0; iev < 1000; ++iev ) {
    vector<Energy> m(np, 0.0*GeV);
    vector<LorentzMomentum> p = SimplePhaseSpace::CMSn(1.0*GeV, m);
    for ( int i = 0; i < np - 1; ++i ) for ( int j = i + 1; j < np; ++j ) {
      Energy2 s = (p[i] + p[j]).m2();
      //      Energy2 kt2 = ktclus2(p[i], p[j]);
      Energy2 kt2 = durham(p[i], p[j]);
      massless << s << '\t' << kt2 << endl;
      for ( int i2 = 0; i2 < np - 1; ++i2 )
	for ( int j2 = i2 + 1; j2 < np; ++j2 ) {
	  if ( i == i2 && j == j2 ) continue;
	  LorentzMomentum p1 = p[i] + p[j];
	  LorentzMomentum p2 = p[i2] + p[j2];
	  Energy2 sm = (p1 + p2).m2()
	    //	    - p1.m2() - p2.m2()
	    ;
	  massive << sm << '\t'
	    //		  << ktclus2(p1, p2) << '\t'
		  << durham(p1, p2) << '\t'
		  << p1.m() << '\t'
		  << p2.m() << endl;
	}
    }
  }
*/
/*
  ofstream massive("massive.dat");
  ofstream massless("massless.dat");
  UseRandom rnd(&(eventGenerator()->random()));
  LorentzMomentum inc(0.0*GeV, 0.0*GeV, 0.5*GeV, 0.5*GeV);
  int np = 7;
  for ( int iev = 0; iev < 1000; ++iev ) {
    vector<Energy> m(np, 0.0*GeV);
    vector<LorentzMomentum> p = SimplePhaseSpace::CMSn(1.0*GeV, m);
    for ( int i = 0; i < np; ++i ) {
      Energy2 t = -(inc - p[i]).m2();
      massless << t << '\t' << p[i].perp2() << endl;
      for ( int j = i + 1; j < np; ++j ) {
	LorentzMomentum pm = p[i] + p[j];
	Energy2 t = -(inc - pm).m2();
	massive << t - pm.m2() << '\t' << pm.perp2() << endl;
      }
    }
  }
*/      
}

void TestDecays::testDecay(tEGPtr eg, long id, int N) {
  
  for ( int i = 0; i < N; ++i ) {
    EventPtr event = new_ptr(Event(PPair()));
    StepPtr firstStep = event->newStep();
    // Create an empty step.

    PPtr p = eg->getParticle(id);
    p->set3Momentum(Momentum3(0.0*GeV, 0.0*GeV, 100.0*GeV));
    firstStep->addParticle(p);

    event = eg->generateEvent(*event);

    if ( i < 1000 ) eg->log() << *event;

  }
}

void TestDecays::testGG(tEGPtr eg, Energy W, int N) {
  
  for ( int i = 0; i < N; ++i ) {
    EventPtr event = new_ptr(Event(PPair()));
    StepPtr firstStep = event->newStep();
    // Create an empty step.

    PPtr p1 = eg->getParticle(ParticleID::g);
    p1->set3Momentum(Momentum3(0.5*W, 0.5*W, 0.0*GeV));
    firstStep->addParticle(p1);
    PPtr p2 = eg->getParticle(ParticleID::g);
    p2->set3Momentum(Momentum3(-0.5*W, -0.5*W, 0.0*GeV));
    firstStep->addParticle(p2);
    p1->colourNeighbour(p2);
    p2->colourNeighbour(p1);

    event = eg->generateEvent(*event);

    if ( i < 10 ) eg->log() << *event;

  }
}

