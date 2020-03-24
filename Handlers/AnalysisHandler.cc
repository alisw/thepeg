// -*- C++ -*-
//
// AnalysisHandler.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AnalysisHandler class.
//

#include "AnalysisHandler.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/UtilityBase.h"

using namespace ThePEG;

bool AnalysisHandler::checkHistogramFactory(bool warn) const {
  if ( generator()->histogramFactory() ) return true;
  if ( warn ) generator()->logWarning(
    NoHistFactory() << "No histogram factory was assigned to the "
    << "EventGenerator, hence no histograms will be produced by "
    << name() << "." << Exception::warning);
  return false;
}


FactoryBase & AnalysisHandler::histogramFactory() {
  return *(generator()->histogramFactory());
}

const FactoryBase & AnalysisHandler::histogramFactory() const {
  return *(generator()->histogramFactory());
}

void AnalysisHandler::normalize(tH1DPtr h, CrossSection unit) const {
  histogramFactory().normalizeToXSec(h, unit);
}

void AnalysisHandler::unitNormalize(tH1DPtr h) const {
  histogramFactory().normalizeToUnity(h);
}

IBPtr AnalysisHandler::clone() const {
  return new_ptr(*this);
}

IBPtr AnalysisHandler::fullclone() const {
  return new_ptr(*this);
}

void AnalysisHandler::analyze(tEventPtr event, long, int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  tcEventPtr cevent = event;
  LorentzRotation r = transform(cevent);
  tPVector particles;
  event->selectFinalState(back_inserter(particles));
  Utilities::transform(particles, r);
  analyze(particles, event->weight());
  for ( int i = 0, N = theSlaves.size(); i < N; ++i )
    theSlaves[i]->analyze(particles, event->weight());
  r.invert();
  Utilities::transform(particles, r);
}

LorentzRotation AnalysisHandler::transform(tEventPtr) const {
  return LorentzRotation();
}

LorentzRotation AnalysisHandler::transform(tcEventPtr) const {
  return LorentzRotation();
}

void AnalysisHandler::analyze(const tPVector & particles) {
  for ( int i = 0, N = particles.size(); i < N; ++i ) analyze(particles[i]);
}

void AnalysisHandler::analyze(const tPVector & particles, double weight) {
  analyze(particles);
  for ( int i = 0, N = particles.size(); i < N; ++i )
    analyze(particles[i], weight);
}

void AnalysisHandler::analyze(tPPtr) {}

void AnalysisHandler::analyze(tPPtr, double) {}

void AnalysisHandler::persistentOutput(PersistentOStream & os) const {
  os << theSlaves;
}

void AnalysisHandler::persistentInput(PersistentIStream & is, int) {
  is >> theSlaves;
}

ClassDescription<AnalysisHandler>
AnalysisHandler::initAnalysisHandler;

void AnalysisHandler::Init() {

  static ClassDocumentation<AnalysisHandler> documentation
    ("The ThePEG::AnalysisHandler class is the base class of all "
     "analysis handlers.");

  static RefVector<AnalysisHandler,AnalysisHandler> interfaceSlaves
    ("Slaves",
     "ThePEG::AnalysisHandler objects to be called for the same extracted "
     "particles as this one.",
     &AnalysisHandler::theSlaves, 0, true, false, true, false);

  interfaceSlaves.rank(10);

}
