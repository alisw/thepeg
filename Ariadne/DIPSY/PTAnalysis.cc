// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PTAnalysis class.
//

#include "PTAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "DipoleEventHandler.h"
#include "EventFiller.h"
#include "DiffractiveEventFiller.h"
#include "ThePEG/Utilities/Current.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

using namespace DIPSY;

PTAnalysis::PTAnalysis() {}

PTAnalysis::~PTAnalysis() {}

void PTAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  if ( !event->primaryCollision() ) return;
  const StepVector & steps = event->primaryCollision()->steps();
  double maxpt = 0.0;
  for ( int is = 0, Ns = min(steps.size(), 2); is < Ns ; ++is ) {
    tH1DPtr h = is > 0? hptG1: hptg1;
    tParticleVector
      gluons(steps[is]->particles().begin(), steps[is]->particles().end());
    for ( int i = 0, N = gluons.size(); i < N; ++i )
      if ( abs(gluons[i]->eta()) < 1.0 ) {
	double pt = gluons[i]->momentum().perp()/GeV;
	h->fill(pt, event->weight());
	if ( pt > maxpt ) maxpt = pt;
      }
  }
  hptm1->fill(maxpt, event->weight());
}

LorentzRotation PTAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void PTAnalysis::analyze(const tPVector & particles, double weight) {
  AnalysisHandler::analyze(particles, weight);
  // Calls analyze() for each particle.
  sumw += weight;
}

void PTAnalysis::analyze(tPPtr p, double weight) {
  if ( !p->data().iCharge() ) return;
  double eta = p->eta();
  double pt = p->momentum().perp()/GeV;
  if ( pt > 5.0 ) heta5->fill(eta, weight);
  if ( eta < -7.0 )
    return;
  else if ( eta < -5.0 )
    hpt75->fill(pt, weight);
  else if ( eta < -3.0 )
    hpt53->fill(pt, weight);
  else if ( eta < -1.0 )
    hpt31->fill(pt, weight);
  else if ( eta < 1.0 )
    hpt11->fill(pt, weight);
  else if ( eta < 3.0 )
    hpt13->fill(pt, weight);
  else if ( eta < 5.0 )
    hpt35->fill(pt, weight);
  else if ( eta < 7.0 )
    hpt57->fill(pt, weight);
  //  sumw += weight;
}

IBPtr PTAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr PTAnalysis::fullclone() const {
  return new_ptr(*this);
}


void PTAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  histogramFactory().initrun();
  histogramFactory().registerClient(this);
  hpt75 = histogramFactory().createHistogram1D
    ("pt75",40,0.0,20.0);
  hpt53 = histogramFactory().createHistogram1D
    ("pt53",40,0.0,20.0);
  hpt31 = histogramFactory().createHistogram1D
    ("pt31",40,0.0,20.0);
  hpt11 = histogramFactory().createHistogram1D
    ("pt11",100,0.0,100.0);
  hptm1 = histogramFactory().createHistogram1D
    ("ptm1",100,0.0,100.0);
  hptg1 = histogramFactory().createHistogram1D
    ("ptg1",100,0.0,100.0);
  hptG1 = histogramFactory().createHistogram1D
    ("ptG1",100,0.0,100.0);
  hpt13 = histogramFactory().createHistogram1D
    ("pt13",40,0.0,20.0);
  hpt35 = histogramFactory().createHistogram1D
    ("pt35",40,0.0,20.0);
  hpt57 = histogramFactory().createHistogram1D
    ("pt57",40,0.0,20.0);
  heta5 = histogramFactory().createHistogram1D
    ("eta5",40,-10.0,10.0);
  sumw = 0.0;
}


void PTAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  if ( sumw <= 0.0 ) return;
  double sc = 2.0/sumw;
  hpt75->scale(sc);
  hpt53->scale(sc);
  hpt31->scale(sc);
  hpt11->scale(sc/2.0);
  hptm1->scale(sc/2.0);
  hptg1->scale(sc/2.0);
  hptG1->scale(sc/2.0);
  hpt13->scale(sc);
  hpt35->scale(sc);
  hpt57->scale(sc);
  heta5->scale(sc);
  


}

void PTAnalysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void PTAnalysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<PTAnalysis,AnalysisHandler>
  describeDIPSYPTAnalysis("DIPSY::PTAnalysis", "PTAnalysis.so");

void PTAnalysis::Init() {

  static ClassDocumentation<PTAnalysis> documentation
    ("There is no documentation for the PTAnalysis class");

}

