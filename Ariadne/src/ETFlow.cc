// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ETFlow class.
//

#include "ETFlow.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/MaxCmp.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/PDT/StandardMatchers.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

using namespace ThePEG;

ETFlow::ETFlow() {}

ETFlow::~ETFlow() {}

void ETFlow::analyze(const tPVector & particles, double weight) {
  AnalysisHandler::analyze(particles);

  // First find the scattered lepton. We assume neutral current and
  // will take the corresponding electron with the largest pt in the
  // lab frame.

  if ( particles.empty() ) return;
  tCollPtr coll = particles[0]->birthStep()->collision();

  tPPtr lepton = coll->incoming().first;
  tPPtr proton = coll->incoming().second;
  if ( !LeptonMatcher::Check(lepton->data()) ) swap(lepton, proton);

  MaxCmp<Energy2, tPPtr> sel;
  for ( int i = 0, N = particles.size(); i < N; ++i )
    if ( particles[i]->id() == lepton->id() )
      sel(particles[i]->momentum().perp2(), particles[i]);
  LorentzMomentum l = sel.index()->momentum();
  LorentzMomentum P = proton->momentum();

  Energy sumfwd = ZERO;
  for ( int i = 0, N = particles.size(); i < N; ++i ) {
    if ( particles[i] == sel.index() ) continue;
    double theta =
      particles[i]->momentum().angle(P)*180.0/Constants::pi;
    if ( theta > 4.4 && theta < 15.0 ) sumfwd += particles[i]->momentum().e();
  }
  if ( sumfwd < 0.5*GeV ) return;

  double theta = l.angle(P)*180.0/Constants::pi;
  if ( theta < 157.0 || theta > 172.5 || l.e() < 14.0*GeV ) return;

  LorentzMomentum k = lepton->momentum() - l;
  Energy2 Q2 = -k.m2();
  Energy2 W2 = (k + P).m2();
  if ( W2 < 3000.0*GeV2 ) return;
  double x = Q2/(Q2 + W2);

  if ( x < 0.001 )
    low += weight;
  else
    hiw += weight;

  LorentzRotation R = Utilities::getBoostToCM(make_pair(k, P));

  for ( int i = 0, N = particles.size(); i < N; ++i ) {
    if ( particles[i] == sel.index() ) continue;
    LorentzMomentum p = R*particles[i]->momentum();
    if ( x < 0.001 )
      lox->fill(p.eta(), p.et()*weight/GeV);
    else
      hix->fill(p.eta(), p.et()*weight/GeV);
  }

}

IBPtr ETFlow::clone() const {
  return new_ptr(*this);
}

IBPtr ETFlow::fullclone() const {
  return new_ptr(*this);
}

void ETFlow::dofinish() {
  AnalysisHandler::dofinish();
  lox->scale(10.0/low);
  hix->scale(10.0/hiw);
}

void ETFlow::doinitrun() {
  AnalysisHandler::doinitrun();
  if ( !checkHistogramFactory(true) ) return;
  histogramFactory().registerClient(this);
  histogramFactory().mkdirs("/ETFLow");
  histogramFactory().cd("/ETFlow");
  lox = histogramFactory().createHistogram1D
    ("Low-x", "Low-x", 200, -10.0, 10.0);
  hix = histogramFactory().createHistogram1D
    ("High-x", "High-x", 200, -10.0, 10.0);
  low = 0.0;
  hiw = 0.0;
}

void ETFlow::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void ETFlow::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<ETFlow,AnalysisHandler>
  describeETFlow("ThePEG::ETFlow", "ETFlow.so");

void ETFlow::Init() {

  static ClassDocumentation<ETFlow> documentation
    ("There is no documentation for the ETFlow class");

}

