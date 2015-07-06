// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Emitter class.
//

#include "Emitter.h"
#include "DipoleState.h"
#include "MECorrBase.h"
#include "ReweightBase.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Emitter.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Current.h"

using namespace Ariadne;

bool Emitter::performEmission() {
  tParPtr emitted = perform();
  if ( finalVeto(this, state(), emitted, lastPT2(), *lastEType()) ) return false;
  return true;
}

double Emitter::
preweight(tcEmiPtr dipole, tcDipoleStatePtr state, long id,
	  const EmissionType & type) const {
  double w = 1.0;
  if ( MECorr() ) w *= MECorr()->preweight(dipole, state, id, type);
  const vector<DipoleRWPtr> & reweighters =
    Current<Ariadne::CascadeHandler>()->reweighters();
  for ( int i = 0, N = reweighters.size(); i < N; ++i )
    w *= reweighters[i]->preweight(dipole, state, id, type);
  return w;
}

double Emitter::
reweight(tcEmiPtr dip, tcDipoleStatePtr state, long id,
	 Energy2 pt2, vector<double> & genVar,
	 const EmissionType & type) const {
  double w = 1.0;
  if ( MECorr() ) w *= MECorr()->reweight(dip, state, id, pt2, genVar, type);
  const vector<DipoleRWPtr> & reweighters =
    Current<Ariadne::CascadeHandler>()->reweighters();
  for ( int i = 0, N = reweighters.size(); i < N; ++i )
    w *= reweighters[i]->reweight(dip, state, id, pt2, genVar, type);
  return w;		    
}

bool Emitter::
finalVeto(tcEmiPtr dipole, tcDipoleStatePtr state, ParPtr parton,
	  Energy2 pt2, const EmissionType & type) const {
  if ( MECorr() && MECorr()->finalVeto(dipole, state, parton, pt2, type) ) return true;
  const vector<DipoleRWPtr> & reweighters =
    Current<Ariadne::CascadeHandler>()->reweighters();
  for ( int i = 0, N = reweighters.size(); i < N; ++i )
    if ( reweighters[i]->finalVeto(dipole, state, parton, pt2, type) ) return true;
  return false;
}

Energy2 Emitter::rndsud(double C, Energy2 pt2max, Energy2 pt2min) const {
  double CN = 1.0/(C*alpha0());
  bool reweight =
    handler()->runningCoupling() == CascadeHandler::externalRunning ||
    handler()->runningCoupling() == CascadeHandler::internalRunning;
  while ( true ) {
    double R = handler()->rnd();
    if ( lambdaQCD2() <= 0.0*GeV2 )
      return CN*log(R) < log(pt2min/pt2max)? -1.0*GeV2: pt2max*pow(R, CN);
    if ( CN*log(R) < log(log(pt2min/lambdaQCD2())/log(pt2max/lambdaQCD2())) )
      return -1.0*GeV2;
    pt2max = lambdaQCD2()*pow(pt2max/lambdaQCD2(),pow(R, CN));
    if ( !reweight ) return pt2max;
    double trueAlpha =
      handler()->runningCoupling() == CascadeHandler::externalRunning?
      handler()->standardModel()->alphaS(pt2max):
      handler()->internalAlphaS()->value(pt2max, handler()->SM());
    double weight = trueAlpha*log(pt2max/lambdaQCD2())/alpha0();
    if ( weight > 1.0 ) handler()->generator()->logWarning(
      WeightException() << "The Ariadne::CascadeHandler '" << handler()->name()
      << "' failed to overestimate the alpha_S specified by the StandardModel "
      << "object. If this hapens too often you should contact the authors."
      << Exception::warning);			   
    if ( handler()->rnd() < weight ) return pt2max;
  }
}

Energy2 Emitter::rndsudEM(double C, Energy2 pt2max, Energy2 pt2min) const {
  bool fixedalpha = handler()->alphaEM0() > 0.0;
  double aEM = fixedalpha? handler()->alphaEM0():
    handler()->standardModel()->alphaEM(pt2max);
  double CN = 1.0/(C*aEM);
  while ( true ) {
    double R = handler()->rnd();
    if ( CN*log(R) < log(pt2min/pt2max) ) return -1.0*GeV2;
    pt2max = pt2max*pow(R, CN);
    double weight = fixedalpha? 1.0:
      handler()->standardModel()->alphaEM(pt2max)/aEM;
    if ( weight > 1.0 ) handler()->generator()->logWarning(
      WeightException() << "The Ariadne::CascadeHandler '" << handler()->name()
      << "' failed to overestimate the alpha_EM specified by the StandardModel "
      << "object. If this hapens too often you should contact the authors."
      << Exception::warning);			   
    if ( handler()->rnd() < weight ) return pt2max;
  }
}

void Emitter::setup(Energy2 pt2max) {
  const StandardModelBase & SM = *handler()->standardModel();
  using namespace ThePEG::Constants;
  switch ( handler()->runningCoupling() ) {
  case CascadeHandler::noRunning:
    theAlpha0 = handler()->alpha0();
    return;
  case CascadeHandler::simpleRunning:
    theLambdaQCD2 = sqr(handler()->lambdaQCD());
    theAlpha0 = 12.0*Constants::pi/(33.0 - 2.0*min(int(SM.Nf(pt2max)),5));
    return;
  default: ;
  }

  Energy2 s1 = sqr(handler()->pTCut());
  double a1 = SM.alphaS(s1);
  Energy2 s2 = pt2max;
  double a2 = SM.alphaS(s2);
  if ( a1 == a2 ) {
    theLambdaQCD2 = -1.0*GeV2;
    theAlpha0 = a1;
  } else {
    theLambdaQCD2 = exp(log(s2/s1)/(1.0 - a1/a2))*s1;
    theAlpha0 = a1*log(s1/lambdaQCD2());
  }
}

bool Emitter::
check(double x1, double x3, double y1, double y2, double y3) const {
  double x2 = 2.0 - x1 - x3;
  if ( x1 < 0.0 || x2 < 0.0 || x3 < 0.0 ) return false;
  double pp1 = sqr(0.5*x1) - y1;
  double pp2 = sqr(0.5*x2) - y2;
  double pp3 = sqr(0.5*x3) - y3;
  if ( pp1 < 0.0 || pp2 < 0.0 || pp3 < 0.0 ||
       2.0*(pp1*pp2 + pp2*pp3 + pp3*pp1) <= sqr(pp1) + sqr(pp2) + sqr(pp3) )
    return false;
  return true;
}

void Emitter::persistentOutput(PersistentOStream & os) const {
  os << ounit(theLastPT2, GeV2) << theAlpha0 << ounit(theLambdaQCD2, GeV2)
     << ounit(theMaxScale, GeV2) << genVar;
}

void Emitter::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theLastPT2, GeV2) >> theAlpha0 >> iunit(theLambdaQCD2, GeV2)
     >> iunit(theMaxScale, GeV2) >> genVar;
  theEType = 0;
}

AbstractClassDescription<Emitter> Emitter::initEmitter;
// Definition of the static class description member.

void Emitter::Init() {}

void Emitter::debugme() const {
  CascadeBase::debugme();
}

bool Emitter::checkIntegrety() {
  return true;
}
