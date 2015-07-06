// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISAnalysis class.
//

#include "DISAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DISAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/RemnantParticle.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/PDT/RemnantDecayer.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/LastXCombInfo.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Utilities/UtilityBase.h"

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

using namespace Ariadne;

using namespace ThePEG;

DISAnalysis::~DISAnalysis() {}

void DISAnalysis::doinit() throw(InitException) {
  AnalysisHandler::doinit();
  checkHistogramFactory(true);
}

void DISAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  if ( !checkHistogramFactory(true) ) return;
  histogramFactory().registerClient(this);
  histogramFactory().mkdirs("/DISAnalysis");
  histogramFactory().cd("/DISAnalysis");
  histY = histogramFactory().createHistogram1D
    ("y", "y", 100, 0.0, 1.0);
  histQ2 = histogramFactory().createHistogram1D
    ("Q2", "Q2", 100, 0.0, 100.0);
  histW = histogramFactory().createHistogram1D
    ("W", "W", 100, 0.0, 100.0);
  histX = histogramFactory().createHistogram1D
    ("logx", "Log10(Bjorken-x)", 100, -4.0, 0.0);
  histYQ = histogramFactory().createHistogram1D
    ("yq", "Quark rapidity", 100, -6.0, 4.0);
  histPTQ = histogramFactory().createHistogram1D
    ("pty", "Quark pT", 40, 0.0, 20.0);
  histYG = histogramFactory().createHistogram1D
    ("yg", "Gluon rapidity", 100, -6.0, 4.0);
  histPTG = histogramFactory().createHistogram1D
    ("ptg", "Gluon pT", 40, 0.0, 20.0);
}

void DISAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  if ( !checkHistogramFactory() ) return;
  normalize(histY);
  normalize(histQ2);
  normalize(histW);
  normalize(histX);
  normalize(histYQ);
  normalize(histPTQ);
  normalize(histYG);
  normalize(histPTG);
}

void DISAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  const PPair & incoming = event->primaryCollision()->incoming();
  tcPPtr proton = incoming.first;
  tcPPtr iniLep = incoming.second;
  if(incoming.second->id() == ParticleID::pplus){
    proton = incoming.second;
    iniLep = incoming.first;
  }

  const SubProcess & sub = *event->primarySubProcess();

  tcPPtr q = sub.outgoing()[0];
  tcPPtr finLep = sub.outgoing()[1];
  if( ! QuarkMatcher::Check(q->data()) ){
    q = sub.outgoing()[1];
    finLep = sub.outgoing()[0];
  }

  LorentzMomentum pIniLep = iniLep->momentum();
  LorentzMomentum pFinLep = finLep->momentum();
  pp = proton->momentum();
  pq = pIniLep - pFinLep;
  Energy2 Q2 = -pq.m2();
  Energy W = (pq + pp).m();
  histQ2->fill( Q2/GeV2 );
  histX->fill(log10( Q2 / (2 * pp.dot(pq)) ));
  histY->fill( pp.dot(pq) / pp.dot(pIniLep) );
  histW->fill( W / GeV );

  AnalysisHandler::analyze(event, ieve, loop, state);
}

LorentzRotation DISAnalysis::transform(tEventPtr event) const {
  return Utilities::getBoostToCM( make_pair( pp, pq ) );
}

void DISAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
}

void DISAnalysis::analyze(tPPtr p) {
    if(QuarkMatcher::Check(p->data())){
      tParticleVector parents = p->parents();
      bool rem = false;
      for ( int j = 0, M = parents.size(); j < M; ++j ) {
        if(dynamic_ptr_cast<tcRemPPtr>(parents[j])){
          rem = true;
        }
      }
      if(!rem){
        histYQ->fill( p->rapidity() );
        histPTQ->fill( p->mt() / GeV );
      }
    }
    else if(p->id() == ParticleID::g){
      histYG->fill( p->rapidity() );
      histPTG->fill( p->mt() / GeV );
    }
}

void DISAnalysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void DISAnalysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<DISAnalysis> DISAnalysis::initDISAnalysis;
// Definition of the static class description member.

void DISAnalysis::Init() {

  static ClassDocumentation<DISAnalysis> documentation
    ("There is no documentation for the DISAnalysis class");

}

