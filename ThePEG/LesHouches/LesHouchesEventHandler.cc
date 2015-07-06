// -*- C++ -*-
//
// LesHouchesEventHandler.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LesHouchesEventHandler class.
//

#include "LesHouchesEventHandler.h"
#include "LesHouchesReader.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Handlers/LuminosityFunction.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/Utilities/LoopGuard.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Utilities/Maths.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

LesHouchesEventHandler::~LesHouchesEventHandler() {}

IBPtr LesHouchesEventHandler::clone() const {
  return new_ptr(*this);
}

IBPtr LesHouchesEventHandler::fullclone() const {
  return new_ptr(*this);
}

void LesHouchesEventHandler::doinit() {
  EventHandler::doinit();
  for ( int i = 0, N = readers().size(); i < N; ++i )
    readers()[i]->init();
}

void LesHouchesEventHandler::initialize() {

  if ( lumiFnPtr() ) Repository::clog()
    << "The LuminosityFunction '" << lumiFnPtr()->name()
    << "' assigned to the LesHouchesEventHandler '" << name()
    << "' will not be active in this run. Instead the incoming "
    << "particles will be determined by the used LesHouchesReader objects.\n"
    << Exception::warning;

  if ( readers().empty() )
    throw LesHouchesInitError()
      << "No readers were defined for the LesHouchesEventHandler '"
      << name() << "'" << Exception::warning;

  // Go through all the readers and collect information about cross
  // sections and processes.
  typedef map<int,tLesHouchesReaderPtr> ProcessMap;
  ProcessMap processes;
  PDPair incoming;
  Energy MaxEA = ZERO;
  Energy MaxEB = ZERO;
  
  for ( int i = 0, N = readers().size(); i < N; ++i ) {
    LesHouchesReader & reader = *readers()[i];
    reader.initialize(*this);

    // Check that the incoming particles are consistent between the
    // readers.
    if ( !incoming.first ) {
      incoming.first = getParticleData(reader.heprup.IDBMUP.first);
      if ( !incoming.first )
	Throw<LesHouchesInitError>() 
	  << "Unknown beam PID " << reader.heprup.IDBMUP.first
	  << ". Have you created a matching BeamParticle object?"
	  << Exception::runerror;
    }
    if ( !incoming.second ) {
      incoming.second = getParticleData(reader.heprup.IDBMUP.second);
      if ( !incoming.second )
	Throw<LesHouchesInitError>() 
	  << "Unknown beam PID " << reader.heprup.IDBMUP.first
	  << ". Have you created a matching BeamParticle object?"
	  << Exception::runerror;
    }
    if ( incoming.first->id() != reader.heprup.IDBMUP.first ||
	 incoming.second->id() != reader.heprup.IDBMUP.second )
      Repository::clog()
	<< "The different LesHouchesReader objects in the "
	<< "LesHouchesEventHandler '" << name() << "' have different "
	<< "types of colliding particles." << Exception::warning;
    MaxEA = max(MaxEA,  reader.heprup.EBMUP.first*GeV);
    MaxEB = max(MaxEB,  reader.heprup.EBMUP.second*GeV);

    // Check that the weighting of the events in the different readers
    // is consistent with the ones requested for this event
    // handler. Also collect the sum of the maximum weights.
    if ( reader.negativeWeights() && weightOption() > 0 )
      throw LesHouchesInitError()
	<< "The reader '" << reader.name()
	<< "' contains negatively weighted events, "
	<< "which is not allowed for the LesHouchesEventHandler '"
	<< name() << "'." << Exception::warning;

    // Check that we do not have the same process numbers in different
    // readers.
    for ( int ip = 0; ip < reader.heprup.NPRUP; ++ip ) {
      if ( reader.heprup.LPRUP[ip] ) {
	ProcessMap::iterator pit = processes.find(reader.heprup.LPRUP[ip]);
	if ( pit == processes.end() )
	  processes[reader.heprup.LPRUP[ip]] = readers()[i];
	else if ( warnPNum ) {
	  Throw<LesHouchesPNumException>()
	    << "In the LesHouchesEventHandler '"
	    << name() << "', both the '" << pit->second->name() << "' and '"
	    << reader.name() << "' contains sub-process number " << pit->first
	    << ". This process may be double-counted in this run."
	    << Exception::warning;
	}
      }
    }
    selector().insert(reader.stats.maxXSec(), i);
  }

  stats.maxXSec(selector().sum());
  histStats.maxXSec(selector().sum());

  // Check that we have any cross section at all.
  if ( stats.maxXSec() <= ZERO )
    throw LesHouchesInitError()
      << "The sum of the cross sections of the readers in the "
      << "LesHouchesEventHandler '" << name()
      << "' was zero." << Exception::warning;

  // We now create a LuminosityFunction object to inform others about
  // the energy of the beam.
  theIncoming = incoming;
  lumiFn(new_ptr(LuminosityFunction(MaxEA, MaxEB)));

}

void LesHouchesEventHandler::doinitrun() {
  EventHandler::doinitrun();
  stats.reset();
  histStats.reset();

  for ( int i = 0, N = readers().size(); i < N; ++i ) {
    readers()[i]->initrun();
  }
}

EventPtr LesHouchesEventHandler::generateEvent() {

  LoopGuard<EventLoopException,LesHouchesEventHandler>
    loopGuard(*this, maxLoop());

  while ( true ) {
    loopGuard();

    currentReader(readers()[selector().select(UseRandom::current())]);

    skipEvents();
    currentReader()->reset();

    double weight = currentReader()->getEvent();
    if ( weightOption() == unitweight && weight < 0.0 ) weight = 0.0;

    if ( weightOption() == unitweight || weightOption() == unitnegweight ) {
      CrossSection newmax = selector().reweight(weight);
      if ( newmax > CrossSection() )
	increaseMaxXSec(newmax);
    }

    select(weight/currentReader()->preweight);
    histStats.select(weight);

    if ( !weighted() ) {
      if ( weightOption() == unitweight  || weightOption() == unitnegweight ) {
	if ( !rndbool(abs(weight)) ) continue;
	weight = Math::sign(1.0, weight);
      }
      else if ( weight == 0.0 ) continue;
    } else if ( weight == 0.0 ) continue;

    accept();

    // Divide by the bias introduced by the preweights in the reader.
    weight /= currentReader()->preweight;

    try {
      theLastXComb = currentReader()->getXComb();

      currentEvent(new_ptr(Event(lastParticles(), this, generator()->runName(),
				 generator()->currentEventNumber(), weight)));
      currentEvent()->optionalWeights() = currentReader()->optionalEventWeights();
      performCollision();

      if ( !currentCollision() ) throw Veto();

      return currentEvent();
    }
    catch (Veto) {
      reject(weight);
    }
    catch (Stop) {
    }
    catch (Exception &) {
      reject(weight);
      throw;
    }
  }
}

void LesHouchesEventHandler::skipEvents() {

  // Don't do this for readers which seem to generate events on the fly.
  if ( currentReader()->active() || currentReader()->NEvents() <= 0 ) return;

  // Estimate the fration of the total events available from
  // currentReader() which will be requested.
  double frac = currentReader()->stats.maxXSec()/stats.maxXSec();
  if ( stats.accepted() > 0 )
    frac *= double(stats.attempts())/double(stats.accepted());
  else
    frac *= double(stats.attempts() + 1);
  double xscan = generator()->N()*frac/currentReader()->NEvents();

  // Estimate the number of times we need to go through the events for
  // the currentReader(), and how many events on average we need to
  // skip for each attempted event to go through the file an integer
  // number of times.
  double nscan = ceil(xscan);
  double meanskip = nscan/xscan - 1.0;

  // Skip an average numer of steps with a Poissonian distribution.
  currentReader()->
    skip(UseRandom::rndPoisson(meanskip)%currentReader()->NEvents());
}

void LesHouchesEventHandler::select(double weight) {
  stats.select(weight);
  currentReader()->select(weight);
}

tCollPtr LesHouchesEventHandler::performCollision() {
  lastExtractor()->select(lastXCombPtr());
  if ( CKKWHandler() ) CKKWHandler()->setXComb(lastXCombPtr());
  currentCollision(new_ptr(Collision(lastParticles(), currentEvent(), this)));
  if ( currentEvent() ) currentEvent()->addCollision(currentCollision());
  currentStep(new_ptr(Step(currentCollision(), this)));
  currentCollision()->addStep(currentStep());
  currentStep()->addSubProcess(currentReader()->getSubProcess());
  lastExtractor()->constructRemnants(lastXCombPtr()->partonBinInstances(),
				     subProcess(), currentStep());

  if ( !currentReader()->cuts().passCuts(*currentCollision()) ) throw Veto();

  initGroups();
  if ( ThePEG_DEBUG_ITEM(1) ) {
    if ( currentEvent() )    
      generator()->logfile() << *currentEvent();
    else
      generator()->logfile() << *currentCollision();
  }
  return continueCollision();
}

EventPtr LesHouchesEventHandler::continueEvent() {

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

void LesHouchesEventHandler::dofinish() {
  EventHandler::dofinish();
  if ( selector().compensating() ) generator()->log()
    << "Warning: The run was ended while the LesHouchesEventHandler '"
    << name() << "' was still trying to compensate for weights larger than 1. "
    << "The cross section estimates may therefore be statistically "
    << "inaccurate." << endl;
}

void LesHouchesEventHandler::statistics(ostream & os) const {
  if ( statLevel() == 0 ) return;

  string line = "======================================="
    "=======================================\n";

  if ( stats.accepted() <= 0 ) {
    os << line << "No events generated by event handler '" << name() << "'."
       << endl;
      return;
  }

  os << line << "Statistics for Les Houches event handler \'" << name() << "\':\n"
     << "                                       "
     << "generated    number of    Cross-section\n"
     << "                                       "
     << "   events     attempts             (nb)\n";

  os << line << "Total:" << setw(42) << stats.accepted() << setw(13)
     << stats.attempts() << setw(17)
     << ouniterr(stats.xSec(), stats.xSecErr(), nanobarn) << endl
     << line;

  if ( statLevel() == 1 ) return;

  if ( statLevel() == 2 ) {

    os << "Per Les Houches Reader breakdown:\n";
  
    for ( int i = 0, N = readers().size(); i < N; ++i ) {
      LesHouchesReader & reader = *readers()[i];
      string n = reader.name();
      n.resize(37, ' ');
      os << n << setw(11) << reader.stats.accepted() << setw(13)
	 << reader.stats.attempts() << setw(17)
	 << ouniterr(reader.stats.xSec(), reader.stats.xSecErr(), nanobarn)
	 << endl;
    }
    os << line;
  } else {

    os << "Per Les Houches Reader (and process #) breakdown:\n";
  
    for ( int i = 0, N = readers().size(); i < N; ++i ) {
      LesHouchesReader & reader = *readers()[i];
      string n = reader.name() + " (all)";
      n.resize(37, ' ');
      os << n << setw(11) << reader.stats.accepted() << setw(13)
	 << reader.stats.attempts() << setw(17)
	 << ouniterr(reader.stats.xSec(), reader.stats.xSecErr(), nanobarn)
	 << endl;
      CrossSection xsectot = reader.stats.xSec();
      if ( xsectot != ZERO ) xsectot /= reader.stats.sumWeights();
      typedef LesHouchesReader::StatMap::const_iterator const_iterator;
      for ( const_iterator i = reader.statmap.begin();
	    i != reader.statmap.end(); ++i ) {
	ostringstream ss;
	ss << reader.name() << " (" << i->first << ")";
	string n = ss.str();
	n.resize(37, ' ');
	os << n << setw(11) << i->second.accepted() << setw(13)
	   << i->second.attempts() << setw(17)
	   << ouniterr(i->second.sumWeights()*xsectot,
		       sqrt(i->second.sumWeights2())*xsectot, nanobarn) << endl;
      }
      os << line;
    }
  }

  string warn = "Warning: Result may be statistically incorrect since\n"
    " the following LesHouchesReaders were oversampled:\n";
  for ( int i = 0, N = readers().size(); i < N; ++i ) {
    LesHouchesReader & reader = *readers()[i];
    if ( reader.NEvents() > 0 && reader.stats.attempts() > reader.NEvents() ) {
      os << warn;
      warn = "";
      os << "'" << reader.name() << "' (by a factor "
	 << double(reader.stats.attempts())/double(reader.NEvents())
	 << ")" << endl;
    }
  }
}


void LesHouchesEventHandler::increaseMaxXSec(CrossSection maxxsec) {
  stats.maxXSec(selector().sum());
  histStats.maxXSec(selector().sum());
  currentReader()->increaseMaxXSec(maxxsec);
}

void LesHouchesEventHandler::accept() {
  stats.accept();
  histStats.accept();
  currentReader()->accept();
}

void LesHouchesEventHandler::reject(double w) {
  stats.reject(w);
  histStats.reject(w);
  currentReader()->reject(w);
}

CrossSection LesHouchesEventHandler::histogramScale() const {
  return histStats.xSec()/histStats.sumWeights();
}

CrossSection LesHouchesEventHandler::integratedXSec() const {
  return histStats.xSec();
}

CrossSection LesHouchesEventHandler::integratedXSecErr() const {
  return histStats.xSecErr();
}

void LesHouchesEventHandler::persistentOutput(PersistentOStream & os) const {
  os << stats << histStats << theReaders << theSelector
     << oenum(theWeightOption) << theUnitTolerance << theCurrentReader << warnPNum;
}

void LesHouchesEventHandler::persistentInput(PersistentIStream & is, int) {
  is >> stats >> histStats >> theReaders >> theSelector
     >> ienum(theWeightOption) >> theUnitTolerance >> theCurrentReader >> warnPNum;
}

ClassDescription<LesHouchesEventHandler>
LesHouchesEventHandler::initLesHouchesEventHandler;
// Definition of the static class description member.

void LesHouchesEventHandler::setUnitTolerance(double x) {
  theUnitTolerance = x;
  selector().tolerance(unitTolerance());
}

void LesHouchesEventHandler::Init() {

  static ClassDocumentation<LesHouchesEventHandler> documentation
    ("This is the main class administrating the selection of hard "
     "subprocesses from a set of ThePEG::LesHouchesReader objects.");


  static RefVector<LesHouchesEventHandler,LesHouchesReader>
    interfaceLesHouchesReaders
    ("LesHouchesReaders",
     "Objects capable of reading events from an event file or an "
     "external matrix element generator.",
     &LesHouchesEventHandler::theReaders, -1, false, false, true, false, false);

  static Switch<LesHouchesEventHandler,WeightOpt> interfaceWeightOption
    ("WeightOption",
     "The different ways to weight events in the Les Houches event handler. "
     "Whether weighted or not and whether or not negative weights are allowed.",
     &LesHouchesEventHandler::theWeightOption, unitweight, true, false);
  static SwitchOption interfaceWeightOptionUnitWeight
    (interfaceWeightOption,
     "UnitWeight",
     "All events have unit weight.",
     unitweight);
  static SwitchOption interfaceWeightOptionNegUnitWeight
    (interfaceWeightOption,
     "NegUnitWeight",
     "All events have weight +1 or maybe -1.",
     unitnegweight);
  static SwitchOption interfaceWeightOptionVarWeight
    (interfaceWeightOption,
     "VarWeight",
     "Events may have varying but positive weights.",
     varweight);
  static SwitchOption interfaceWeightOptionVarNegWeight
    (interfaceWeightOption,
     "VarNegWeight",
     "Events may have varying weights, both positive and negative.",
     varnegweight);

  static Switch<LesHouchesEventHandler,bool> interfaceWarnPNum
    ("WarnPNum",
     "Warn if the same process number is used in more than one "
     "LesHouchesReader.",
     &LesHouchesEventHandler::warnPNum, true, true, false);
  static SwitchOption interfaceWarnPNumWarning
    (interfaceWarnPNum,
     "Warning",
     "Give a warning message.",
     true);
  static SwitchOption interfaceWarnPNumNoWarning
    (interfaceWarnPNum,
     "NoWarning",
     "Don't give a warning message.",
     false);

  static Parameter<LesHouchesEventHandler,double> interfaceUnitTolerance
    ("UnitTolerance",
     "If the <interface>WeightOption</interface> is set to unit weight, do not start compensating unless the a weight is found to be this much larger than unity.",
     &LesHouchesEventHandler::theUnitTolerance, 1.0e-6, 0.0, 0,
     true, false, Interface::lowerlim,
     &LesHouchesEventHandler::setUnitTolerance,
     (double(LesHouchesEventHandler::*)()const)(0),
     (double(LesHouchesEventHandler::*)()const)(0),
     (double(LesHouchesEventHandler::*)()const)(0),
     (double(LesHouchesEventHandler::*)()const)(0));

  interfaceLesHouchesReaders.rank(10);
  interfaceWeightOption.rank(9);

}

