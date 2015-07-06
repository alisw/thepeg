// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GapAnalysis class.
//

#include "GapAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Utilities/UtilityBase.h"

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

GapAnalysis::GapAnalysis(): ok(true) {}

GapAnalysis::~GapAnalysis() {}

void GapAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  if ( !ok ) return;
  if ( event->incoming().second->id() != ParticleID::pplus ) {
    ok = false;
    Throw<Exception>() << "DIPSY::GapAnalysis did not find a proton as secons incoming particle "
		       << "The Analysis will be switched off." << Exception::warning;
    return;
  }

  Lorentz5Momentum pp = event->incoming().second->momentum();
  Energy Ep = 820.0*GeV;
  Lorentz5Momentum ppnew = Lorentz5Momentum(ZERO, ZERO, -sqrt(sqr(Ep) - pp.m2()), Ep, pp.m());
  LorentzRotation R = Utilities::transformToMomentum(pp, ppnew);
  tPVector particles;
  event->selectFinalState(back_inserter(particles));

  double etamax = -1000.0;
  LorentzMomentum sumx;
  for ( int i = 0, N = particles.size(); i < N; ++i ) {
    if ( particles[i]->momentum().perp() < 30.0*MeV ) continue;
    double eta = -(R*particles[i]->momentum()).eta();
    if ( eta < 4.3 && eta > etamax ) etamax = eta;
    if ( eta < 4.3 ) sumx += particles[i]->momentum();
  }
  if ( etamax == -1000.0 )
    Throw<Exception>() << "DIPSY::GapAnalysis found an empty detector." << Exception::warning;

  histEtamax->fill(etamax, event->weight());
  histMX->fill(sumx.m()/GeV, event->weight());

}

LorentzRotation GapAnalysis::transform(tcEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void GapAnalysis::analyze(const tPVector & particles, double weight) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void GapAnalysis::analyze(tPPtr, double weight) {}

void GapAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  histogramFactory().normalizeToUnity(histEtamax);
  histogramFactory().normalizeToUnity(histMX);
}

void GapAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  histogramFactory().registerClient(this); // Initialize histograms.
  histogramFactory().mkdirs("/HZ93093"); // Put histograms in specal directory.
  histEtamax = histogramFactory().createHistogram1D("/HZ93093/etamax", 20, -5.0, 5.0);
  histMX = histogramFactory().createHistogram1D("/HZ93093/MX", 25, 0.0, 100.0);
}


IBPtr GapAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr GapAnalysis::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void GapAnalysis::persistentOutput(PersistentOStream & os) const {
  os << oenum(ok); // Add all member variable which should be written persistently here.
}

void GapAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> ienum(ok);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<GapAnalysis,AnalysisHandler>
  describeDIPSYGapAnalysis("DIPSY::GapAnalysis", "GapAnalysis.so");

void GapAnalysis::Init() {

  static ClassDocumentation<GapAnalysis> documentation
    ("There is no documentation for the GapAnalysis class");

}

