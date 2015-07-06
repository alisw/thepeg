// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HIAnalysis class.
//

#include "HIAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Utilities/Throw.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

using namespace DIPSY;

HIAnalysis::HIAnalysis(): tuple(0) {}

HIAnalysis::~HIAnalysis() {}

void HIAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {

  LorentzPoint p1 = event->incoming().first->vertex();
  LorentzPoint p2 = event->incoming().second->vertex();

  double b = (p1 - p2).perp()/femtometer;

  *tuple << event->weight() << '\t' << b << '\t';

  AnalysisHandler::analyze(event, ieve, loop, state);
}

LorentzRotation HIAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void HIAnalysis::analyze(const tPVector & particles, double weight) {
  //  AnalysisHandler::analyze(particles, weight);
  // Calls analyze() for each particle.
  int Nch = 0;
  double sumETFwd = 0.0;
  LorentzPoint sample;
  for ( int i = 0, N = particles.size(); i < N; ++i ) {
    const Particle & pa = *particles[i];
    const ParticleData & pd = pa.data();
    const Lorentz5Momentum & p = pa.momentum();
    try {
      double etabs = abs(p.eta());
      if ( etabs > 3.2 && etabs <= 4.9 ) sumETFwd += p.et()/GeV;
      if ( etabs < 0.5 && pd.iCharge() ) {
	sample = pa.vertex();
	++Nch;
	double pt = p.perp()/GeV;
	ptch->fill(pt, weight);
	if ( abs(pd.id()) == ParticleID::piplus ) ptpi->fill(pt, weight);
	if ( abs(pd.id()) == ParticleID::Kplus ) ptK->fill(pt, weight);
	if ( abs(pd.id()) == ParticleID::pplus ) ptp->fill(pt, weight);
      }
    } catch ( Exception & e ) {
      e.handle();
      Throw<Exception>() << e.message() << Exception::warning;
    }
  }
  *tuple << Nch << '\t' << sumETFwd << endl;
}

IBPtr HIAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr HIAnalysis::fullclone() const {
  return new_ptr(*this);
}


void HIAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  if ( tuple ) delete tuple;
  tuple = new ofstream((generator()->runName() + ".tuple").c_str());
  histogramFactory().initrun();
  histogramFactory().registerClient(this);
  ptch = histogramFactory().createHistogram1D
    ("ptch", 100, 0.0, 10.0);
  ptpi = histogramFactory().createHistogram1D
    ("ptpi", 100, 0.0, 10.0);
  ptK = histogramFactory().createHistogram1D
    ("ptK", 100, 0.0, 10.0);
  ptp = histogramFactory().createHistogram1D
    ("ptp", 100, 0.0, 10.0);
}


void HIAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  if ( tuple ) delete tuple;
  tuple = 0;
  double sumw = generator()->sumWeights();
  if ( sumw > 0.0 ) {
    ptch->scale(10.0/sumw);
    ptpi->scale(10.0/sumw);
    ptK->scale(10.0/sumw);
    ptp->scale(10.0/sumw);
  }
}

void HIAnalysis::persistentOutput(PersistentOStream & os) const {
  //  os <<; 
}

void HIAnalysis::persistentInput(PersistentIStream & is, int) {
  //  is >>;
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<HIAnalysis,AnalysisHandler>
  describeDIPSYHIAnalysis("DIPSY::HIAnalysis", "HIAnalysis.so");

void HIAnalysis::Init() {

  static ClassDocumentation<HIAnalysis> documentation
    ("There is no documentation for the HIAnalysis class");

}

