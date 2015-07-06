// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DYAnalysis class.
//

#include "DYAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/EventGenerator.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DYAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

using namespace ThePEG;

using namespace Ariadne;

DYAnalysis::~DYAnalysis() {}

void DYAnalysis::doinit() throw(InitException) {
  AnalysisHandler::doinit();
  checkHistogramFactory(true);
}

void DYAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  if ( !checkHistogramFactory(true) ) return;
  histogramFactory().registerClient(this);
  histogramFactory().mkdirs("/DYAnalysis");
  histogramFactory().cd("/DYAnalysis");
  histPt = histogramFactory().createHistogram1D
    ("Pt", "Transverse momentum of Drell-Yan pair", 100, 0.0, 400.0);
  histEta = histogramFactory().createHistogram1D
    ("Eta", "Pseudo-rapidity of Drell-Yan pair", 40, -5.0, 5.0);
  histY = histogramFactory().createHistogram1D
    ("Y", "Rapidity of Drell-Yan pair", 40, -5.0, 5.0);
  histMass = histogramFactory().createHistogram1D
    ("Mass", "Invariant mass of Drell-Yan pair", 100, 0.0, 500.0);
  histCPt = histogramFactory().createHistogram1D
    ("CPt", "Transverse momentum of Drell-Yan pair in |y|&lt;1",
     100, 0.0, 400.0);
  histYG = histogramFactory().createHistogram1D
    ("YG", "Rapidity of the emitted gluon", 40, -5.0, 5.0);
  histYQ = histogramFactory().createHistogram1D
    ("YQ", "Rapidity of the emitted quark", 40, -5.0, 5.0);
  histNG = histogramFactory().createHistogram1D
    ("NG", "Number of gluons", 50, 0.0, 50.0);
  histPTG = histogramFactory().createHistogram1D
    ("PTG", "Transverse momentum of gluons", 100, 0.0, 100.0);
  histPTQ = histogramFactory().createHistogram1D
    ("PTQ", "Transverse momentum of quarks", 100, 0.0, 100.0);
}

void DYAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  if ( !checkHistogramFactory() ) return;
  normalize(histPt);
  normalize(histEta);
  normalize(histY);
  normalize(histYG);
  normalize(histYQ);
  normalize(histMass);
  normalize(histCPt);
  normalize(histNG);
  normalize(histPTG);
  normalize(histPTQ);
}

void DYAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
}

LorentzRotation DYAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void DYAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  if ( !checkHistogramFactory() ) return;
  tPPair dypair;
  Energy2 mins = Constants::MaxEnergy2;
  int ng = 0;
  for ( int i = 0, N = particles.size(); i < N; ++i ) {
    if ( particles[i]->id() == ParticleID::g ){
      histYG->fill(particles[i]->momentum().rapidity());
      ng++;
      if(abs(particles[i]->momentum().eta()) < 2.5){
        histPTG->fill(particles[i]->momentum().perp()/GeV);
      }
    }

    if ( QuarkMatcher::Check(particles[i]->data()) ) {
      LorentzMomentum pq = particles[i]->momentum();
      histPTQ->fill(pq.perp()/GeV);
      if(pq.perp() > 20*GeV){
        histYQ->fill(pq.rapidity());
      }
    }

    if ( !LeptonMatcher::Check(particles[i]->data()) ) continue;
    if ( !ChargedMatcher::Check(particles[i]->data()) ) continue;
    for ( int j = i + 1, M = particles.size(); j < M; ++j ) {
      if ( particles[i]->id() + particles[j]->id() != 0 ) continue;
      Energy2 s = ( particles[i]->momentum() + particles[j]->momentum() ).m2();
      if ( s >= mins ) continue;
      dypair = make_pair(particles[i], particles[j]);
      mins = s;
    }
  }
  histNG->fill(ng);
  if ( !dypair.first ) return;
  LorentzMomentum pdy = dypair.first->momentum() + dypair.second->momentum();
  double y = pdy.rapidity();
  Energy pt = pdy.perp();
  histPt->fill(pt/GeV);
  histEta->fill(pdy.eta());
  histY->fill(y);
  histMass->fill(sqrt(mins)/GeV);
  if ( abs(y) < 1.0 ) histCPt->fill(pt/GeV);
}

void DYAnalysis::analyze(tPPtr) {}

void DYAnalysis::persistentOutput(PersistentOStream & os) const {}

void DYAnalysis::persistentInput(PersistentIStream & is, int) {}

ClassDescription<DYAnalysis> DYAnalysis::initDYAnalysis;
// Definition of the static class description member.

void DYAnalysis::Init() {

  static ClassDocumentation<DYAnalysis> documentation
    ("The DYAnalysis class collects a number of histograms related to "
     "Drell-Yan production.");

}

