// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TestJunctions class.
//

#include "TestJunctions.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TestJunctions.tcc"
#endif


using namespace ThePEG;

TestJunctions::~TestJunctions() {}

NoPIOClassDescription<TestJunctions> TestJunctions::initTestJunctions;
// Definition of the static class description member.

void TestJunctions::Init() {
  breakThePEG();
  if ( N() <= 0 ) N(eventGenerator()->N());
  testSimpleG(eventGenerator(), N()/4); 
  testConnectedG(eventGenerator(), N()/4); 
  testConnected(eventGenerator(), N()/4); 
  testSimple(eventGenerator(), N()/4); 
}

void TestJunctions::testConnected(tEGPtr eg, int N) {
  
  for ( int i = 0; i < N; ++i ) {
    // Create an empty step.
    EventPtr event = new_ptr(Event(PPair()));
    StepPtr firstStep = event->newStep();

    PPtr q1 = eg->getParticle(ParticleID::u);
    q1->set3Momentum(Momentum3(20.0*GeV, ZERO, 100.0*GeV));
    tColinePtr c1 = ColourLine::create(q1);
    PPtr q2 = eg->getParticle(ParticleID::u);
    q2->set3Momentum(Momentum3(-10.0*GeV, ZERO, 50.0*GeV));
    tColinePtr c2 = ColourLine::create(q2);

    PPtr aq1 = eg->getParticle(ParticleID::ubar);
    aq1->set3Momentum(Momentum3(10.0*GeV, ZERO, -100.0*GeV));
    tColinePtr ac1 = ColourLine::createAnti(aq1);
    PPtr aq2 = eg->getParticle(ParticleID::ubar);
    aq2->set3Momentum(Momentum3(-20.0*GeV, ZERO, -50.0*GeV));
    tColinePtr ac2 = ColourLine::createAnti(aq2);

    ColourLine::create(c1, c2, ac1, ac2);

    firstStep->addParticle(q1);
    firstStep->addParticle(q2);
    firstStep->addParticle(aq1);
    firstStep->addParticle(aq2);
    event = eg->generateEvent(*event);

    if ( i < 10 ) eg->log() << *event;

  }
}

void TestJunctions::testConnectedG(tEGPtr eg, int N) {
  
  for ( int i = 0; i < N; ++i ) {
    // Create an empty step.
    EventPtr event = new_ptr(Event(PPair()));
    StepPtr firstStep = event->newStep();

    PPtr q1 = eg->getParticle(ParticleID::u);
    q1->set3Momentum(Momentum3(20.0*GeV, ZERO, 80.0*GeV));
    PPtr g3 = eg->getParticle(ParticleID::g);
    g3->set3Momentum(Momentum3(ZERO, 1.0*GeV, 10.0*GeV));
    PPtr g4 = eg->getParticle(ParticleID::g);
    g4->set3Momentum(Momentum3(ZERO, -1.0*GeV, 10.0*GeV));
    ColourLine::create(tPPtr(q1), tPPtr(g3));
    ColourLine::create(tPPtr(g3), tPPtr(g4));
    tColinePtr c1 = ColourLine::create(g4);
    PPtr q2 = eg->getParticle(ParticleID::u);
    q2->set3Momentum(Momentum3(-10.0*GeV, ZERO, 50.0*GeV));
    tColinePtr c2 = ColourLine::create(q2);

    PPtr aq1 = eg->getParticle(ParticleID::ubar);
    aq1->set3Momentum(Momentum3(10.0*GeV, ZERO, -100.0*GeV));
    tColinePtr ac1 = ColourLine::createAnti(aq1);
    PPtr aq2 = eg->getParticle(ParticleID::ubar);
    aq2->set3Momentum(Momentum3(-20.0*GeV, ZERO, -50.0*GeV));
    tColinePtr ac2 = ColourLine::createAnti(aq2);

    PPtr g1 = eg->getParticle(ParticleID::g);
    g1->set3Momentum(Momentum3(2.0*GeV, ZERO, ZERO));
    PPtr g2 = eg->getParticle(ParticleID::g);
    g2->set3Momentum(Momentum3(-2.0*GeV, ZERO, ZERO));
    tColinePtr cj1 = ColourLine::create(g1);
    tColinePtr cj2 = ColourLine::createAnti(g2);
    ColourLine::create(tPPtr(g2), tPPtr(g1));
    cj1->setSourceNeighbours(c1, c2);
    cj2->setSinkNeighbours(ac1, ac2);
    
    firstStep->addParticle(q1);
    firstStep->addParticle(q2);
    firstStep->addParticle(aq1);
    firstStep->addParticle(aq2);
    firstStep->addParticle(g1);
    firstStep->addParticle(g2);
    firstStep->addParticle(g3);
    firstStep->addParticle(g4);

    event = eg->generateEvent(*event);

    if ( i < 10 ) eg->log() << *event;

  }
}

void TestJunctions::testSimple(tEGPtr eg, int N) {
  
  for ( int i = 0; i < N; ++i ) {
    // Create an empty step.
    EventPtr event = new_ptr(Event(PPair()));
    StepPtr firstStep = event->newStep();

    PPtr q1 = eg->getParticle(ParticleID::u);
    q1->set3Momentum(Momentum3(ZERO, ZERO, 100.0*GeV));
    PPtr q2 = eg->getParticle(ParticleID::u);
    q2->set3Momentum(Momentum3(-20.0*GeV, ZERO, -50.0*GeV));
    PPtr q3 = eg->getParticle(ParticleID::u);
    q3->set3Momentum(Momentum3(20.0*GeV, ZERO, -50.0*GeV));
    ColourLine::create(q1)->setSourceNeighbours(ColourLine::create(q2),
						ColourLine::create(q3));
    firstStep->addParticle(q1);
    firstStep->addParticle(q2);
    firstStep->addParticle(q3);

    event = eg->generateEvent(*event);

    if ( i < 10 ) eg->log() << *event;

  }
}

void TestJunctions::testSimpleG(tEGPtr eg, int N) {
  
  for ( int i = 0; i < N; ++i ) {
    // Create an empty step.
    EventPtr event = new_ptr(Event(PPair()));
    StepPtr firstStep = event->newStep();

    PPtr q1 = eg->getParticle(ParticleID::u);
    q1->set3Momentum(Momentum3(ZERO, ZERO, 100.0*GeV));
    PPtr q2 = eg->getParticle(ParticleID::u);
    q2->set3Momentum(Momentum3(-20.0*GeV, ZERO, -50.0*GeV));
    PPtr q3 = eg->getParticle(ParticleID::u);
    q3->set3Momentum(Momentum3(20.0*GeV, ZERO, -50.0*GeV));
    PPtr g1 = eg->getParticle(ParticleID::g);
    g1->set3Momentum(Momentum3(2.0*GeV, ZERO, ZERO));
    PPtr g2 = eg->getParticle(ParticleID::g);
    g2->set3Momentum(Momentum3(-2.0*GeV, ZERO, ZERO));
    ColourLine::create(tPPtr(q3), tPPtr(g1));
    ColourLine::create(tPPtr(g1), tPPtr(g2));
    ColourLine::create(q1)->setSourceNeighbours(ColourLine::create(q2),
						ColourLine::create(g2));
    firstStep->addParticle(q1);
    firstStep->addParticle(q2);
    firstStep->addParticle(g1);
    firstStep->addParticle(g2);
    firstStep->addParticle(q3);

    event = eg->generateEvent(*event);

    if ( i < 10 ) eg->log() << *event;

  }
}

