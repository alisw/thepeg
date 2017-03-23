// -*- C++ -*-
//
// StandardEventHandler.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardEventHandler class.
//

#include "StandardEventHandler.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Handlers/StdXCombGroup.h"
#include "ThePEG/Handlers/SubProcessHandler.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/LoopGuard.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/MatrixElement/MEGroup.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/LuminosityFunction.h"
#include "ThePEG/Handlers/SamplerBase.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Config/algorithm.h"
#include <iomanip>
#include <sstream>

using namespace ThePEG;

StandardEventHandler::StandardEventHandler()
  : EventHandler(false), collisionCuts(true), theLumiDim(0) {
  setupGroups();
}

StandardEventHandler::~StandardEventHandler() {}

void StandardEventHandler::reject(double w) {
  tStdXCombPtr last = dynamic_ptr_cast<tStdXCombPtr>(lastXCombPtr());
  if ( !last ) return;
  last->reject(w);
  xSecStats.reject(w);
}

void StandardEventHandler::reweight(double factor) const {
  tStdXCombPtr last = dynamic_ptr_cast<tStdXCombPtr>(lastXCombPtr());
  if ( !currentEvent() || !last )
    return;
  double weight = currentEvent()->weight();
  last->reweight(weight,factor*weight);
  xSecStats.reweight(weight,factor*weight);
  currentEvent()->weight(factor*weight);
  for ( map<string,double>::iterator w = 
	  currentEvent()->optionalWeights().begin();
	w != currentEvent()->optionalWeights().end(); ++w )
    w->second *= factor;
}

void StandardEventHandler::doupdate() {
  EventHandler::doupdate();
  bool redo = touched();
  UpdateChecker::check(theIncomingA, redo);
  UpdateChecker::check(theIncomingB, redo);
  for_each(subProcesses(), UpdateChecker(redo));
  if ( !redo ) return;
  theIncoming.first = theIncomingA;
  theIncoming.second = theIncomingB;
  for ( SubHandlerList::iterator sit = subProcesses().begin();
	sit != subProcesses().end(); ++sit )
    if ( !(**sit).pExtractor()->canHandle(incoming()) )
      throw StandardEventHandlerUpdateException()
	<<  "Cannot use the parton extractor '" << (**sit).pExtractor()->name()
	<< "' in the SubProcessHandler '" << (**sit).name() << "' in the "
	<< "StandardEventHandler '" << name() << "' since it cannot handle "
	<< "the requested types of incoming particles ("
	<< theIncomingA->name() << "," << theIncomingB->name() << ").";
}

void StandardEventHandler::doinit() {
  EventHandler::doinit();
  if ( !lumiFnPtr() ) throw StandardEventHandlerUpdateException()
    << "The StandardEventHandler '" << name() << "' does not have any "
    << "LuminosityFunction object assigned to it, which it needs to be "
    << "able to generate events." << Exception::warning;
}

IBPtr StandardEventHandler::clone() const {
  return new_ptr(*this);
}

IBPtr StandardEventHandler::fullclone() const {
  return new_ptr(*this);
}

void StandardEventHandler::
addME(Energy maxEnergy, tSubHdlPtr sub, tPExtrPtr extractor, tCutsPtr cuts,
      tCascHdlPtr ckkw, tMEPtr me, const PBPair & pBins,
      const PartonPairVec& allPBins) {
  typedef MEBase::DiagramVector DiagramVector;
  typedef map<string,DiagramVector> DiagramMap;
  cPDPair pin(pBins.first->parton(), pBins.second->parton());
  DiagramVector diag = me->diagrams();
  DiagramMap tdiag;
  DiagramMap tmdiag;
  for ( int i = 0, N = diag.size(); i < N; ++i ) {
    cPDPair din(diag[i]->partons()[0], diag[i]->partons()[1]);
    if (!me->noMirror())
      if ( din.first->id() < din.second->id() ) swap(din.first, din.second); 
    if ( din == pin ) tdiag[diag[i]->getTag()].push_back(diag[i]);
    if (!me->noMirror())
      if ( din.first == pin.second && din.second == pin.first )
	tmdiag[diag[i]->getTag()].push_back(diag[i]);
  }

  if ( tdiag.empty() ) tdiag = tmdiag;
  for ( DiagramMap::iterator dit = tdiag.begin(); dit != tdiag.end(); ++dit ) {
    cPDPair din(dit->second.back()->partons()[0],
		dit->second.back()->partons()[1]);
    // check
    assert(me->noMirror() ? din == pin : true);
    StdXCombPtr xcomb = me->makeXComb(maxEnergy, incoming(), this, sub, extractor,
				      ckkw, pBins, cuts, dit->second, din != pin,
				      allPBins);
    if ( xcomb->checkInit() ) xCombs().push_back(xcomb);
    else generator()->logWarning(
      StandardEventHandlerInitError() << "The matrix element '"
      << xcomb->matrixElement()->name() 
      << "' cannot generate the diagram '"
      << dit->first 
      << "' when used together with the parton extractor '"
      << xcomb->pExtractor()->name()
      << "'.\nThe corresponding diagram is switched off, " 
      << "check the collision energy and/or the cuts."
      << Exception::warning);
  }
}

void StandardEventHandler::initGroups() {
  tStdXCombPtr lastXC = dynamic_ptr_cast<tStdXCombPtr>(lastXCombPtr());
  if ( lastXC ) optGroups = lastXC->subProcessHandler()->groups();
  EventHandler::initGroups();
}

tCollPtr StandardEventHandler::performCollision() {
  tStdXCombPtr lastXC = dynamic_ptr_cast<tStdXCombPtr>(lastXCombPtr());
  if ( CKKWHandler() ) CKKWHandler()->setXComb(lastXCombPtr());
  lastExtractor()->select(lastXC);
  currentCollision(new_ptr(Collision(lastParticles(), currentEvent(), this)));
  if ( currentEvent() ) currentEvent()->addCollision(currentCollision());
  currentStep(new_ptr(Step(currentCollision())));
  currentCollision()->addStep(currentStep());

  currentStep()->addSubProcess(lastXC->construct());

  lastExtractor()->construct(lastXC->partonBinInstances(), currentStep());
  if ( collisionCuts )
    if ( !lastCuts().passCuts(*currentCollision()) ) throw Veto();
  initGroups();
  if ( ThePEG_DEBUG_ITEM(1) ) {
    if ( currentEvent() )    
      generator()->logfile() << *currentEvent();
    else
      generator()->logfile() << *currentCollision();
  }
  return continueCollision();
}

void StandardEventHandler::setScale(Energy2 scale) {
  lastXCombPtr()->lastScale(scale);
}

void StandardEventHandler::initialize() {

  theLumiDim = lumiFn().nDim(incoming());
  Energy maxEnergy = lumiFn().maximumCMEnergy();

  xCombs().clear();

  cuts()->initialize(sqr(maxEnergy), lumiFn().Y());

  for ( SubHandlerList::const_iterator sit = subProcesses().begin();
	sit != subProcesses().end(); ++sit ) {
    CutsPtr kincuts = (**sit).cuts()? (**sit).cuts(): cuts();
    if ( (**sit).cuts() ) kincuts->initialize(sqr(maxEnergy), lumiFn().Y());
    PExtrPtr pextract = (**sit).pExtractor();

    tCascHdlPtr ckkw = (**sit).CKKWHandler();
    if ( !ckkw ) ckkw = CKKWHandler();

    PartonPairVec vpc = pextract->getPartons(maxEnergy, incoming(), *kincuts);

    // The last parton bin pair was in fact the bins corresponding to
    // the incoming particles, so we remove them, but save them to
    // keep them referenced.
    PBPair orig = vpc.back();
    vpc.pop_back();
    for ( PartonPairVec::iterator ppit = vpc.begin();
	  ppit != vpc.end(); ++ppit )
      for ( MEVector::const_iterator meit = (**sit).MEs().begin();
	    meit != (**sit).MEs().end(); ++meit ) {
	addME(maxEnergy, *sit, pextract, kincuts, ckkw, *meit, *ppit,vpc);
      }
  }

  theMaxDims.clear();
  for ( int i = 0, N = xCombs().size(); i < N; ++i )
    theMaxDims.push_back(xCombs()[i]->nDim());

  sampler()->setEventHandler(this);
  sampler()->initialize();
}

CrossSection StandardEventHandler::
dSigDR(const pair<double,double> ll, Energy2 maxS,
       int ibin, int nr, const double * r) {
  PPair inc = make_pair(incoming().first->produceParticle(),
			incoming().second->produceParticle());
  SimplePhaseSpace::CMS(inc, maxS);
  xCombs()[ibin]->prepare(inc);
  return xCombs()[ibin]->dSigDR(ll, nr, r);
}

tStdXCombPtr StandardEventHandler::select(int bin, double & weight) {
  tStdXCombPtr lastXC = xCombs()[bin];
  // clean up the old XComb object before switching to a new one
  if ( theLastXComb && theLastXComb != lastXC ) theLastXComb->clean();
  theLastXComb = lastXC;
  weight /= lastXC->matrixElement()->preWeight();
  lastXC->select(weight);
  xSecStats.select(weight);
  lastXC->accept();
  xSecStats.accept();
  return lastXC;
}

int StandardEventHandler::nBins() const {
  return xCombs().size();
}
 
void StandardEventHandler::statistics(ostream & os) const {
  if ( statLevel() == 0 ) return;
  map<cPDPair, XSecStat> partonMap;
  map<MEPtr, XSecStat> meMap;
  map<PExtrPtr, XSecStat> extractMap;
  XSecStat tot(sampler()->maxXSec());

  for ( int i = 0, N = xCombs().size(); i < N; ++i ) {
    const StandardXComb & x = *xCombs()[i];
    if ( partonMap.find(x.partons()) == partonMap.end() )
      partonMap[x.partons()] = XSecStat(sampler()->maxXSec());
    partonMap[x.partons()] += x.stats();
    if ( meMap.find(x.matrixElement()) == meMap.end() )
      meMap[x.matrixElement()] = XSecStat(sampler()->maxXSec());
    meMap[x.matrixElement()] += x.stats();
    if ( extractMap.find(x.pExtractor()) == extractMap.end() )
      extractMap[x.pExtractor()] = XSecStat(sampler()->maxXSec());
    extractMap[x.pExtractor()] += x.stats();
    tot += x.stats();
  }

  string line = string(78, '=') + "\n";

  if ( tot.accepted() <= 0 ) {
    os << line << "No events generated by event handler '" << name() << "'."
       << endl;
      return;
  }

  os << line << "Statistics for event handler \'" << name() << "\':\n"
     << "                                       "
     << "generated    number of    Cross-section\n"
     << "                                       "
     << "   events     attempts             (nb)\n";

  os << line << "Total (from   weighted events): including vetoed events" << setw(23)
     << ouniterr(sampler()->integratedXSec(), 
		 sampler()->integratedXSecErr(), nanobarn)
     << endl;
  os << line << "Total (from "
     << (weighted() ? "  weighted" : "unweighted") << " events):" 
     << setw(17) << tot.accepted() << setw(13)
     << tot.attempts() << setw(17)
     << ouniterr(tot.xSec(sampler()->attempts()),tot.xSecErr(sampler()->attempts()) , nanobarn)
     << endl << line;

  if ( statLevel() == 1 ) return;

  os << "Per matrix element breakdown:\n";
  for ( map<MEPtr, XSecStat>::iterator i = meMap.begin();
	i != meMap.end(); ++i ) {
    string n = i->first->name();
    n.resize(37, ' ');
    os << n << setw(11) << i->second.accepted() << setw(13)
       << i->second.attempts() << setw(17)
       << ouniterr(i->second.xSec(sampler()->attempts()), i->second.xSecErr(sampler()->attempts()), nanobarn)
       << endl;
  }
  os << line;

  if ( statLevel() == 2 ) return;

  os << "Per parton extractor breakdown:\n";
  for ( map<PExtrPtr, XSecStat>::iterator i = extractMap.begin();
	i != extractMap.end(); ++i ) {
    string n = i->first->name();
    n.resize(37, ' ');
    os << n << setw(11) << i->second.accepted() << setw(13)
       << i->second.attempts() << setw(17)
       << ouniterr(i->second.xSec(sampler()->attempts()), i->second.xSecErr(sampler()->attempts()), nanobarn)
       << endl;
  }
  os << line;

  os << "Per incoming partons breakdown:\n";
  for ( map<cPDPair, XSecStat>::iterator i = partonMap.begin();
	i != partonMap.end(); ++i ) {
    string n = i->first.first->PDGName() + " " + i->first.second->PDGName();
    n.resize(37, ' ');
    os << n << setw(11) << i->second.accepted() << setw(13)
       << i->second.attempts() << setw(17)
       << ouniterr(i->second.xSec(sampler()->attempts()), i->second.xSecErr(sampler()->attempts()), nanobarn)
       << endl;
  }
  os << line;

  if ( statLevel() == 3 ) return;

  os << "Detailed breakdown:\n";

  for ( int i = 0, N = xCombs().size(); i < N; ++i ) {
    const StandardXComb & x = *xCombs()[i];
    XSecStat xstat(sampler()->maxXSec());
    xstat += x.stats();
    os << "(" << x.pExtractor()->name() << ") "
       << x.partons().first->PDGName() << " "
       << x.partons().second->PDGName()
       << " (" << x.matrixElement()->name() << " "
       << x.lastDiagram()->getTag() << ") " << endl
       << setw(48) << xstat.accepted() << setw(13) << xstat.attempts()
       << setw(17)
       << ouniterr(xstat.xSec(sampler()->attempts()), xstat.xSecErr(sampler()->attempts()), nanobarn) << endl;
  }

  os << line;

}

CrossSection StandardEventHandler::histogramScale() const {
  xSecStats.maxXSec(sampler()->maxXSec());
  return xSecStats.xSec(sampler()->attempts())/xSecStats.sumWeights();
}

CrossSection StandardEventHandler::integratedXSec() const {
  xSecStats.maxXSec(sampler()->maxXSec());
  return xSecStats.xSec(sampler()->attempts());
}

CrossSection StandardEventHandler::integratedXSecErr() const {
  xSecStats.maxXSec(sampler()->maxXSec());
  return xSecStats.xSecErr(sampler()->attempts());
}

CrossSection StandardEventHandler::integratedXSecNoReweight() const {
  xSecStats.maxXSec(sampler()->maxXSec());
  return xSecStats.xSecNoReweight(sampler()->attempts());
}

CrossSection StandardEventHandler::integratedXSecErrNoReweight() const {
  xSecStats.maxXSec(sampler()->maxXSec());
  return xSecStats.xSecErrNoReweight(sampler()->attempts());
}

void StandardEventHandler::doinitrun() {
  EventHandler::doinitrun();
  for ( SubHandlerList::iterator sit = subProcesses().begin();
	sit != subProcesses().end(); ++sit )
    (**sit).initrun();
  sampler()->initrun();
  for ( int i = 0, N = xCombs().size(); i < N; ++i )
    xCombs()[i]->reset();
  xSecStats.reset();
}

CrossSection StandardEventHandler::dSigDR(const vector<double> & r) {
  double jac = 1.0;
  pair<double,double> ll = lumiFn().generateLL(&r[0], jac);
  Energy2 maxS = sqr(lumiFn().maximumCMEnergy())/exp(ll.first + ll.second);
  int bin = sampler()->lastBin();
  CrossSection x = jac*lumiFn().value(incoming(), ll.first, ll.second)
    *dSigDR(ll, maxS, bin, nDim(bin) - lumiDim(), &r[lumiDim()]);
  return x;
}

EventPtr StandardEventHandler::generateEvent() {

  LoopGuard<EventLoopException,StandardEventHandler>
    loopGuard(*this, maxLoop());

  while (1) {
    loopGuard();

    EventHandler::clean();

    double weight = sampler()->generate();

    tStdXCombPtr lastXC = 
      select(sampler()->lastBin(), weight);

    try {

      lumiFn().select(lastXC);

      currentEventBoost() = lumiFn().getBoost();

      currentEvent(new_ptr(Event(lastParticles(), this, generator()->runName(),
				 generator()->currentEventNumber(), weight)));

      performCollision();
      if ( !currentCollision() ) throw Veto();

      currentEvent()->transform(currentEventBoost());

      return currentEvent();

    }
    catch (Veto) {
      reject(currentEvent()->weight());
    }
    catch (Stop) {
      break;
    }
    catch (Exception &) {
      reject(currentEvent()->weight());
      throw;
    }
  }
  return currentEvent();
}

EventPtr StandardEventHandler::continueEvent() {
  if ( !generator() ) throw StandardEventHandlerInitError()
    << "The event handler '" << name() << "' had not been isolated "
    << "from the setup phase before it was used." << Exception::maybeabort;
  try {
    continueCollision();
  }
  catch (Veto) {
    reject(currentEvent()->weight());
  }
  catch (Stop) {
  }
  catch (Exception &) {
    reject(currentEvent()->weight());
    throw;
  }
  return currentEvent(); 
}

void StandardEventHandler::select(tXCombPtr newXComb) {
  EventHandler::select(newXComb);
  lastExtractor()->select(newXComb);
}

void StandardEventHandler::clean() {
  if ( theLastXComb ) theLastXComb->clean();
  for (size_t i=0; i < theXCombs.size(); ++i )
    if ( theXCombs[i] ) theXCombs[i]->clean();
  EventHandler::clean();
}

void StandardEventHandler::dofinish() {
  clean();
  EventHandler::dofinish();
}

ClassDescription<StandardEventHandler>
StandardEventHandler::initStandardEventHandler;	// 

void StandardEventHandler::setIncomingA(PDPtr p) {
  theIncomingA = p;
  theIncoming.first = p;
}

void StandardEventHandler::setIncomingB(PDPtr p) {
  theIncomingB = p;
  theIncoming.second = p;
}

void StandardEventHandler::Init() {

  static ClassDocumentation<StandardEventHandler> documentation
    ("This is the standard event handler to generate hard sub-processes "
     "within ThePEG. It must specify a pair of incoming particle beams "
     "in <interface>BeamA</interface> and <interface>BeamB</interface> "
     "and a suiteable <interface>LuminosityFunction</interface>. In "
     "addition at least one object describing the sub-processes to be "
     "generated must be specified in "
     "<interface>SubProcessHandlers</interface>.");

  static Reference<StandardEventHandler,ParticleData> interfaceIncomingA
    ("BeamA",
     "The type of particles in first beam",
     &StandardEventHandler::theIncomingA, false, false, true, false,
     &StandardEventHandler::setIncomingA);

  static Reference<StandardEventHandler,ParticleData> interfaceIncomingB
    ("BeamB",
     "The type of particles in second beam",
     &StandardEventHandler::theIncomingB, false, false, true, false,
     &StandardEventHandler::setIncomingB);

  static RefVector<StandardEventHandler,SubProcessHandler> interfaceSubhandlers
    ("SubProcessHandlers",
     "The list of sub-process handlers used in this StandardEventHandler. ",
     &StandardEventHandler::theSubProcesses, 0, false, false, true, false);

  static Reference<StandardEventHandler,Cuts> interfaceCuts
    ("Cuts",
     "Common kinematical cuts for this StandardEventHandler. These cuts "
     "may be overidden in individual sub-process handlers.",
     &StandardEventHandler::theCuts, false, false, true, false);

  static Switch<StandardEventHandler,bool> interfaceCollisionCuts
    ("CollisionCuts",
     "Switch on or off cuts on collision objects",
     &StandardEventHandler::collisionCuts, true, false, false);
  static SwitchOption interfaceCollisionCutsOn
    (interfaceCollisionCuts,
     "On",
     "Switch on cuts on collision objects",
     true);
  static SwitchOption interfaceCollisionCutsOff
    (interfaceCollisionCuts,
     "Off",
     "Switch off cuts on collision cuts",
     false);

  static Reference<StandardEventHandler,SamplerBase> interfaceSampler
    ("Sampler",
     "The phase space sampler responsible for generating phase space"
     "points according to the cross section given by this event handler",
     &StandardEventHandler::theSampler, false, false, true, true);

  interfaceSubhandlers.rank(11);
  interfaceIncomingA.rank(3);
  interfaceIncomingB.rank(2);

}

void StandardEventHandler::persistentOutput(PersistentOStream & os) const {
  os << theIncomingA << theIncomingB << theSubProcesses << theCuts << collisionCuts
     << theXCombs << theMaxDims << theSampler << theLumiDim << xSecStats;
}

void StandardEventHandler::persistentInput(PersistentIStream & is, int) {
  is >> theIncomingA >> theIncomingB >> theSubProcesses >> theCuts >> collisionCuts
     >> theXCombs >> theMaxDims >> theSampler >> theLumiDim >> xSecStats;
}

