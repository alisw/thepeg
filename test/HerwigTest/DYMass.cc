// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DYMass class.
//

#include "DYMass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/EventGenerator.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DYMass.tcc"
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
using namespace ThePEG;

DYMass::~DYMass() {}

void DYMass::doinitrun() {
  AnalysisHandler::doinitrun();
  histogramFactory().registerClient(this);
  histogramFactory().mkdirs("/DYMass");
  hist =
    histogramFactory().createHistogram1D("/DYMass/dymass", 100, 0.0, 5000.0);
  hist->setTitle("d\\sigma/dm (pb/GeV)");
  vector<double> edges;
  edges.push_back(1100.0);
  edges.push_back(1200.0);
  edges.push_back(1300.0);
  edges.push_back(1500.0);
  edges.push_back(2000.0);
  edges.push_back(3000.0);
  vhist =
    histogramFactory().createHistogram1D("/DYMass/vdymass", "", edges);
  vhist->setTitle("d\\sigma/dm (pb/GeV) [peak]");
  histogramFactory().createDataSet
    ("/DYMass/dymassdata", "d\\sigma/dm (pb/GeV) [peak]", 2)
      << 1250 << 0 << 0 << 0.00011 << 0.00003  << 0.00003
      << 1400 << 0 << 0 << 0.00008 << 0.00002  << 0.00002
      << 1750 << 0 << 0 << 0.00004 << 0.00001  << 0.00001
      << 2500 << 0 << 0 << 0.00001 << 0.000005 << 0.000005;
}

void DYMass::dofinish() {
  AnalysisHandler::dofinish();
  tH1DPtr hist2 =
    histogramFactory().createHistogram1D("/dummy2", 100, 0.0, 5000.0);
  tH1DPtr hist3 =
    histogramFactory().histogramFactory().add("/dummy3", *hist, *hist2);
  tH1DPtr hist4 =
    histogramFactory().histogramFactory().divide("/dummy4", *hist, *hist3);
  hist4 =
    histogramFactory().histogramFactory().divide("/dummy5", *hist, *hist2);
  hist4 =
    histogramFactory().histogramFactory().multiply("/dummy6", *hist, *hist3);
  normalize(hist);
  normalize(vhist);
}

void DYMass::analyze(tEventPtr event, long ieve, int loop, int state) {
  // AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  tSubProPtr sub = event->primaryCollision()->primarySubProcess();
  if ( sub->outgoing().size() != 2 ) return;
  if ( sub->outgoing()[0]->id() + sub->outgoing()[1]->id() ) return;
  if ( abs(sub->outgoing()[0]->id()) != ParticleID::eminus ) return;
  LorentzMomentum p = sub->outgoing()[0]->momentum() +
    sub->outgoing()[1]->momentum();
  hist->fill(p.m()/GeV, event->weight());
  vhist->fill(p.m()/GeV, event->weight());
}

LorentzRotation DYMass::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void DYMass::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void DYMass::analyze(tPPtr) {}

void DYMass::persistentOutput(PersistentOStream & os) const {}

void DYMass::persistentInput(PersistentIStream & is, int) {}

ClassDescription<DYMass> DYMass::initDYMass;
// Definition of the static class description member.

void DYMass::Init() {

  static ClassDocumentation<DYMass> documentation
    ("There is no documentation for the DYMass class");

}

