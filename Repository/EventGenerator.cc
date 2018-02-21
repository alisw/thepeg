// -*- C++ -*-
//
// EventGenerator.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EventGenerator class.
//

#include "EventGenerator.h"
#include "EventGenerator.xh"
#include "ThePEG/Handlers/EventHandler.h"
#include "Repository.h"
#include "ThePEG/Utilities/HoldFlag.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/DebugItem.h"
#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/Strategy.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Analysis/FactoryBase.h"
#include "ThePEG/Handlers/EventManipulator.h"
#include "ThePEG/Handlers/LuminosityFunction.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Handlers/SubProcessHandler.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/Handlers/HadronizationHandler.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/algorithm.h"
#include "ThePEG/Utilities/DynamicLoader.h"
#include <cstdlib>
#include "ThePEG/Repository/Main.h"
#include <csignal>

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "EventGenerator.tcc"
#endif

using namespace ThePEG;

namespace {
  volatile sig_atomic_t THEPEG_SIGNAL_STATE = 0;
}

// signal handler function
// very restricted in what it is allowed do
// without causing undefined behaviour
extern "C" {
  void thepegSignalHandler(int id) {
    THEPEG_SIGNAL_STATE=id;
    signal(id,SIG_DFL);
  }
}

void EventGenerator::checkSignalState() {
  if (THEPEG_SIGNAL_STATE) {
    log() << "Caught signal " << THEPEG_SIGNAL_STATE << ". Exiting ..." << std::endl;
    finalize();
    exit(0);
  }
}

EventGenerator::EventGenerator()
  : thePath("."), theNumberOfEvents(1000), theQuickSize(7000),
    preinitializing(false), ieve(0), weightSum(0.0),
    theDebugLevel(0), logNonDefault(-1), printEvent(0), dumpPeriod(0),
    keepAllDumps(false),
    debugEvent(0), maxWarnings(10), maxErrors(10), theCurrentRandom(0),
    theCurrentGenerator(0), useStdout(false), theIntermediateOutput(false) {}

EventGenerator::EventGenerator(const EventGenerator & eg)
  : Interfaced(eg), theDefaultObjects(eg.theDefaultObjects),
    theLocalParticles(eg.theLocalParticles),
    theStandardModel(eg.theStandardModel),
    theStrategy(eg.theStrategy), theRandom(eg.theRandom),
    theEventHandler(eg.theEventHandler),
    theAnalysisHandlers(eg.theAnalysisHandlers),
    theHistogramFactory(eg.theHistogramFactory),
    theEventManipulator(eg.theEventManipulator),
    thePath(eg.thePath), theRunName(eg.theRunName),
    theNumberOfEvents(eg.theNumberOfEvents), theObjects(eg.theObjects),
    theObjectMap(eg.theObjectMap),
    theParticles(eg.theParticles), theQuickParticles(eg.theQuickParticles),
    theQuickSize(eg.theQuickSize), preinitializing(false),
    theMatchers(eg.theMatchers),
    usedObjects(eg.usedObjects), ieve(eg.ieve), weightSum(eg.weightSum),
    theDebugLevel(eg.theDebugLevel), logNonDefault(eg.logNonDefault),
    printEvent(eg.printEvent), dumpPeriod(eg.dumpPeriod),
    keepAllDumps(eg.keepAllDumps),
    debugEvent(eg.debugEvent),
    maxWarnings(eg.maxWarnings), maxErrors(eg.maxErrors), theCurrentRandom(0),
    theCurrentGenerator(0),
    theCurrentEventHandler(eg.theCurrentEventHandler),
    theCurrentStepHandler(eg.theCurrentStepHandler),
    useStdout(eg.useStdout),
    theIntermediateOutput(eg.theIntermediateOutput) {}

EventGenerator::~EventGenerator() {
  if ( theCurrentRandom ) delete theCurrentRandom;
  if ( theCurrentGenerator ) delete theCurrentGenerator;
}

IBPtr EventGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr EventGenerator::fullclone() const {
  return new_ptr(*this);
}

tcEventPtr EventGenerator::currentEvent() const {
  return eventHandler()->currentEvent();
}

CrossSection EventGenerator::histogramScale() const {
  return eventHandler()->histogramScale();
}

CrossSection EventGenerator::integratedXSec() const {
  return eventHandler()->integratedXSec();
}

CrossSection EventGenerator::integratedXSecErr() const {
  return eventHandler()->integratedXSecErr();
}

void EventGenerator::setSeed(long seed) {
  random().setSeed(seed);
  ostringstream s;
  s << seed;
  const InterfaceBase * ifb = BaseRepository::FindInterface(theRandom, "Seed");
  ifb->exec(*theRandom, "set", s.str());
}

void
EventGenerator::setup(string newRunName,
		      ObjectSet & newObjects,
		      ParticleMap & newParticles,
		      MatcherSet & newMatchers) {
  HoldFlag<int> debug(Debug::level, Debug::isset? Debug::level: theDebugLevel);
  theRunName = newRunName;
  theObjects.swap(newObjects);
  theParticles.swap(newParticles);
  theMatchers.swap(newMatchers);
  theObjectMap.clear();
  for ( ObjectSet::const_iterator it = objects().begin();
	it != objects().end(); ++it ) theObjectMap[(**it).fullName()] = *it;
  UseRandom currentRandom(theRandom);
  CurrentGenerator currentGenerator(this);

  // Force update of all objects and then reset.
  touch();
  for_each(theObjects, mem_fun(&InterfacedBase::touch));
  update();
  for_each(theObjects, mem_fun(&InterfacedBase::update));
  clear();
  BaseRepository::clearAll(theObjects);

  init();

}

IBPtr EventGenerator::getPointer(string name) const {
  ObjectMap::const_iterator it = objectMap().find(name);
  if ( it == objectMap().end() ) return IBPtr();
  else return it->second;
}

void EventGenerator::openOutputFiles() {
  if ( !useStdout ) {
    logfile().open((filename() + ".log").c_str());
    theOutFileName = filename() + ".out";
    outfile().open(theOutFileName.c_str());
    outfile().close();
    theOutStream.str("");
  }
  out() << Repository::banner() << endl;
  log() << Repository::banner() << endl;
}

void EventGenerator::closeOutputFiles() {
  flushOutputFile();
  if ( !useStdout ) logfile().close();
}

void EventGenerator::flushOutputFile() {
  if ( !useStdout ) {
    outfile().open(theOutFileName.c_str(), ios::out|ios::app);
    outfile() << theOutStream.str();
    outfile().close();
  } else
    BaseRepository::cout() << theOutStream.str();
  theOutStream.str("");
}

void EventGenerator::doinit() {

  HoldFlag<int> debug(Debug::level, Debug::isset? Debug::level: theDebugLevel);

  // First initialize base class and random number generator.
  Interfaced::doinit();
  random().init();

  // Make random generator and this available in standard static
  // classes.
  UseRandom useRandom(theRandom);
  CurrentGenerator currentGenerator(this);

  // First initialize all objects which have requested this by
  // implementing a InterfacedBase::preInitialize() function which
  // returns true.
  while ( true ) {
    HoldFlag<bool> hold(preinitializing, true);
    ObjectSet preinits;
    for ( ObjectSet::iterator it = objects().begin();
	  it != objects().end(); ++it )
      if ( (**it).preInitialize() &&
	   (**it).state() == InterfacedBase::uninitialized )
	preinits.insert(*it);
    if ( preinits.empty() ) break;
    for_each(preinits, mem_fun(&InterfacedBase::init));
  }

  // Initialize the quick access to particles.
  theQuickParticles.clear();
  theQuickParticles.resize(2*theQuickSize);
  for ( ParticleMap::const_iterator pit = theParticles.begin();
	pit != theParticles.end(); ++pit )
    if ( abs(pit->second->id()) < theQuickSize )
      theQuickParticles[pit->second->id()+theQuickSize] = pit->second;

  // Then call the init method for all objects. Start with the
  // standard model and the strategy.
  standardModel()->init();
  if ( strategy() ) strategy()->init();
  eventHandler()->init();

  // initialize particles first
  for(ParticleMap::const_iterator pit = particles().begin();
      pit != particles().end(); ++pit) pit->second->init();
  
  for_each(objects(), mem_fun(&InterfacedBase::init));

  // Then initialize the Event Handler calculating initial cross
  // sections and stuff.
  eventHandler()->initialize();
}

void EventGenerator::doinitrun() {

  HoldFlag<int> debug(Debug::level, Debug::isset? Debug::level: theDebugLevel);

  signal(SIGHUP, thepegSignalHandler);
  signal(SIGINT, thepegSignalHandler);
  signal(SIGTERM,thepegSignalHandler);

  currentEventHandler(eventHandler());

  Interfaced::doinitrun();
  random().initrun();

  // Then call the init method for all objects. Start with the
  // standard model and the strategy.
  standardModel()->initrun();
  if ( strategy() ) { 
    strategy()->initrun();
    if ( ! strategy()->versionstring().empty() ) {
      out() << ">> " << strategy()->versionstring() << '\n' << endl;
      log() << ">> " << strategy()->versionstring() << '\n' << endl;
    }
  }
  // initialize particles first
  for(ParticleMap::const_iterator pit = particles().begin();
      pit != particles().end(); ++pit) {
    pit->second->initrun();
  }
  eventHandler()->initrun();

  
  for_each(objects(), mem_fun(&InterfacedBase::initrun));

  if ( logNonDefault > 0 || ( ThePEG_DEBUG_LEVEL && logNonDefault == 0 ) ) {
    vector< pair<IBPtr, const InterfaceBase *> > changed =
      Repository::getNonDefaultInterfaces(objects());
    if ( changed.size() ) {
      log() << string(78, '=') << endl
	    << "The following interfaces have non-default values (default):"
	    << endl << string(78, '-') << endl;
      for ( int i = 0, N = changed.size(); i < N; ++i ) {
	log() << changed[i].first->fullName() << ":"
	      << changed[i].second->name() << " = "
	      << changed[i].second->exec(*changed[i].first, "notdef", "")
	      << endl;
      }
      log() << string(78,'=') << endl;
    }
  }

  weightSum = 0.0;

}

PDPtr EventGenerator::getParticleData(PID id) const {
  long newId = id;
  if ( abs(newId) < theQuickSize && theQuickParticles.size() )
    return theQuickParticles[newId+theQuickSize];
  ParticleMap::const_iterator it = theParticles.find(newId);
  if ( it == theParticles.end() ) return PDPtr();
  return it->second;
}

PPtr EventGenerator::getParticle(PID newId) const {
  tcPDPtr pd = getParticleData(newId);
  if ( !pd ) return PPtr();
  return pd->produceParticle();
}

void EventGenerator::finalize() {
  UseRandom currentRandom(theRandom);
  CurrentGenerator currentGenerator(this);
  finish();
  finally();
}

void EventGenerator::dofinish() {

  HoldFlag<int> debug(Debug::level, Debug::isset? Debug::level: theDebugLevel);

  // first write out statistics from the event handler.
  eventHandler()->statistics(out());

  // Call the finish method for all other objects.
  for_each(objects(), mem_fun(&InterfacedBase::finish));

  if ( theExceptions.empty() ) {
    log() << "No exceptions reported in this run.\n";
  } else {

    log() << "\nThe following exception classes were reported in this run:\n";

    for ( ExceptionMap::iterator it = theExceptions.begin();
	  it != theExceptions.end(); ++it ) {
      string severity;
      switch ( it->first.second ) {
      case Exception::info       : severity="info"; break;
      case Exception::warning    : severity="warning"; break;
      case Exception::setuperror : severity="setuperror"; break;
      case Exception::eventerror : severity="eventerror"; break;
      case Exception::runerror   : severity="runerror"; break;
      case Exception::maybeabort : severity="maybeabort"; break;
      case Exception::abortnow   : severity="abortnow"; break;
      default                    : severity="unknown";
      }
      log() << it->first.first << ' ' << severity
	    << " (" << it->second << " times)\n";
    }
  }

  theExceptions.clear();

  const string & msg = theMiscStream.str();
  if ( ! msg.empty() ) {
    log() << endl 
	  << "Miscellaneous output from modules to the standard output:\n\n"
	  << msg;
    theMiscStream.str("");
  }

  flushOutputFile();

}

void EventGenerator::finally() {

  generateReferences();

  closeOutputFiles();

  if ( theCurrentRandom ) delete theCurrentRandom;
  if ( theCurrentGenerator ) delete theCurrentGenerator;
  theCurrentRandom = 0;
  theCurrentGenerator = 0;

}

void EventGenerator::initialize(bool initOnly) {
  UseRandom currentRandom(theRandom);
  CurrentGenerator currentGenerator(this);
  doInitialize(initOnly);
}

bool EventGenerator::loadMain(string file) {
  initialize();
  UseRandom currentRandom(theRandom);
  CurrentGenerator currentGenerator(this);
  Main::eventGenerator(this);
  bool ok = DynamicLoader::load(file);
  finish();
  finally();
  return ok;
}


void EventGenerator::go(long next, long maxevent, bool tics) {
  UseRandom currentRandom(theRandom);
  CurrentGenerator currentGenerator(this);
  doGo(next, maxevent, tics);
}

EventPtr EventGenerator::shoot() {
  static DebugItem debugfpu("ThePEG::FPU", 1);
  if ( debugfpu ) Debug::unmaskFpuErrors();
  UseRandom currentRandom(theRandom);
  CurrentGenerator currentGenerator(this);
  checkSignalState();
  EventPtr event = doShoot();
  if ( event ) weightSum += event->weight();
  DebugItem::tic();
  return event;
}

EventPtr EventGenerator::doShoot() {
  EventPtr event;
  if ( N() >= 0 && ++ieve > N() ) return event;
  HoldFlag<int> debug(Debug::level, Debug::isset? Debug::level: theDebugLevel);
  do { 
    int state = 0;
    int loop = 1;
    eventHandler()->clearEvent();
    try {
      do {
	// Generate a full event or part of an event
	if ( eventHandler()->empty() ) event = eventHandler()->generateEvent();
	else event = eventHandler()->continueEvent();

	if ( eventHandler()->empty() ) loop = -loop;
	
	// Analyze the possibly uncomplete event
	for ( AnalysisVector::iterator it = analysisHandlers().begin();
	      it != analysisHandlers().end(); ++it )
	  (**it).analyze(event, ieve, loop, state);
	
	// Manipulate the current event, possibly deleting some steps
	// and telling the event handler to redo them.
	if ( manipulator() )
	  state = manipulator()->manipulate(eventHandler(), event);
	
	// If the event was not completed, continue generation and continue.
	loop = abs(loop) + 1;
      } while ( !eventHandler()->empty() );
    }
    catch (Exception & ex) {
      if ( logException(ex, eventHandler()->currentEvent()) ) throw;
    }
    catch (...) {
      event = eventHandler()->currentEvent();
      if ( event )
	log() << *event;
      else
	log() << "An exception occurred before any event object was created!";
      log() << endl;
      dump();
      throw;
    }
    if ( ThePEG_DEBUG_LEVEL ) {
      if ( ( ThePEG_DEBUG_LEVEL == Debug::printEveryEvent ||
	     ieve < printEvent ) && event ) log() << *event;
      if ( debugEvent > 0 && ieve + 1 >= debugEvent )
	Debug::level = Debug::full;
    }
  } while ( !event );

  // If scheduled, dump a clean state between events
  if ( ThePEG_DEBUG_LEVEL && dumpPeriod > 0 && ieve%dumpPeriod == 0 ) {
    eventHandler()->clearEvent();
    eventHandler()->clean();
    dump();
  }

  return event;
}

EventPtr EventGenerator::doGenerateEvent(tEventPtr e) {
  if ( N() >= 0 && ++ieve > N() ) return EventPtr();
  EventPtr event = e;
  try {
    event = eventHandler()->generateEvent(e);
  }
  catch (Exception & ex) {
    if ( logException(ex, eventHandler()->currentEvent()) ) throw;
  }
  catch (...) {
    event = eventHandler()->currentEvent();
    if ( !event ) event = e;
    log() << *event << endl;
    dump();
    throw;
  }
  return event;
}

EventPtr EventGenerator::doGenerateEvent(tStepPtr s) {
  if ( N() >= 0 && ++ieve > N() ) return EventPtr();
  EventPtr event;
  try {
    event = eventHandler()->generateEvent(s);
  }
  catch (Exception & ex) {
    if ( logException(ex, eventHandler()->currentEvent()) ) throw;
  }
  catch (...) {
    event = eventHandler()->currentEvent();
    if ( event ) log() << *event << endl;
    dump();
    throw;
  }
  return event;
}

EventPtr EventGenerator::generateEvent(Event & e) {
  UseRandom currentRandom(theRandom);
  CurrentGenerator currentGenerator(this);
  EventPtr event = doGenerateEvent(tEventPtr(&e));
  if ( event ) weightSum += event->weight();
  return event;
}

EventPtr EventGenerator::generateEvent(Step & s) {
  UseRandom currentRandom(theRandom);
  CurrentGenerator currentGenerator(this);
  EventPtr event = doGenerateEvent(tStepPtr(&s));
  if ( event ) weightSum += event->weight();
  return event;
}

Energy EventGenerator::maximumCMEnergy() const {
  tcEHPtr eh = eventHandler();
  return eh->lumiFnPtr()? eh->lumiFn().maximumCMEnergy(): ZERO;
}

void EventGenerator::doInitialize(bool initOnly) {

  if ( !initOnly )
    openOutputFiles();

  init();

  if ( !initOnly )
    initrun();

  if ( !ThePEG_DEBUG_LEVEL ) Exception::noabort = true;

}

void EventGenerator::doGo(long next, long maxevent, bool tics) {

  if ( maxevent >= 0 ) N(maxevent);

  if ( next >= 0 ) {
    if ( tics ) 
      cerr << "event> " << setw(9) << "init\r" << flush;
    initialize();
    ieve = next-1;
  } else {
    openOutputFiles();
  }

  if ( tics ) tic();
  try {
    while ( shoot() ) {
      if ( tics ) tic();
    }
  }
  catch ( ... ) {
    finish();
    throw;
  }

  finish();

  finally();

}

void EventGenerator::tic(long currev, long totev) const {
  if ( !currev ) currev = ieve;
  if ( !totev ) totev = N();
  long i = currev;
  long n = totev;
  bool skip = currev%(max(totev/100, 1L));
  if ( i > n/2 ) i = n-i;
  while ( skip && i >= 10 && !(i%10) ) i /= 10;
  if ( i == 1 || i == 2 || i == 5 ) skip = false;
  if (!theIntermediateOutput) { //default
    if ( skip ) return;
    cerr << "event> " << setw(8) << currev << " " << setw(8) << totev << "\r";
  }
  else if (theIntermediateOutput) {
    if ( skip && currev%10000!=0) return;
    cerr <<  "event> " << setw(9) << right << currev << "/" << totev 
         << "; xs = " << integratedXSec()/picobarn << " pb +- " 
         << integratedXSecErr()/picobarn << " pb" << endl;
  }
  cerr.flush();
  if ( currev == totev ) cerr << endl;
}
  

void EventGenerator::dump() const {
  if ( dumpPeriod > -1 ) {
    string dumpfile;
    if ( keepAllDumps ) {
      ostringstream number;
      number << ieve;
      dumpfile = filename() + "-" + number.str() + ".dump";
    }
    else
      dumpfile = filename() + ".dump";
    PersistentOStream file(dumpfile, globalLibraries());
    file << tcEGPtr(this);
  }
}

void EventGenerator::use(const Interfaced & i) {
  IBPtr ip = getPtr(i);
  if ( ip ) usedObjects.insert(ip);
}

void EventGenerator::generateReferences() {
  typedef map<string,string> StringMap;
  StringMap references;

  // First get all model descriptions and model references from the
  // used objects. Put them in a map indexed by the description to
  // avoid duplicates.
  for ( ObjectSet::iterator it = usedObjects.begin();
	it != usedObjects.end(); ++it ) {
    if ( *it == strategy() ) continue;
    string desc = Repository::getModelDescription(*it);
    if ( desc.empty() ) continue;
    if ( dynamic_ptr_cast<cEHPtr>(*it) ) desc = "A " + desc;
    else if ( dynamic_ptr_cast<cSMPtr>(*it) ) desc = "B " + desc;
    else if ( dynamic_ptr_cast<cMEPtr>(*it) ) desc = "C " + desc;
    else if ( dynamic_ptr_cast<cCascHdlPtr>(*it) ) desc = "D " + desc;
    else if ( dynamic_ptr_cast<cHadrHdlPtr>(*it) ) desc = "E " + desc;
    else if ( dynamic_ptr_cast<cStepHdlPtr>(*it) ) desc = "F " + desc;
    else if ( dynamic_ptr_cast<cDecayerPtr>(*it) ) desc = "Y " + desc;
    else if ( dynamic_ptr_cast<cAnalysisHdlPtr>(*it) ) desc = "Z " + desc;
    else if ( dynamic_ptr_cast<Ptr<HandlerBase>::const_pointer>(*it) )
      desc = "G " + desc;
    else  desc = "H " + desc;
    references[desc] = Repository::getModelReferences(*it);
  }

  // Now get the main strategy description which should put first and
  // remove it from the map.
  string stratdesc;
  string stratref;
  if ( strategy() ) {
    stratdesc = Repository::getModelDescription(strategy());
    stratref = Repository::getModelReferences(strategy());
    references.erase(stratdesc);
  }

  // Open the file and write out an appendix header
  if ( !useStdout )
    reffile().open((filename() + ".tex").c_str());
  ref() << "\\documentclass{article}\n"
	<< "\\usepackage{graphics}\n"
	<< "\\begin{document}\n"
	<< "\\appendix\n"
	<< "\\section[xxx]{\\textsc{ThePEG} version " << Repository::version()
	<< " \\cite{ThePEG} Run Information}\n"
	<< "Run name: \\textbf{" << runName()
	<< "}:\\\\\n";
  if ( !stratdesc.empty() )
    ref() << "This run was generated using " << stratdesc 
	  << " and the following models:\n";
  else
    ref() << "The following models were used:\n";

  ref() << "\\begin{itemize}\n";

  // Write out all descriptions.
  for ( StringMap::iterator it = references.begin();
	it != references.end(); ++it )
    ref() << "\\item " << it->first.substr(2) << endl;

  // Write out thebibliography header and all references.
  ref() << "\\end{itemize}\n\n"
	    << "\\begin{thebibliography}{99}\n"
	    << "\\bibitem{ThePEG} L.~L\\\"onnblad, "
	    << "Comput.~Phys.~Commun.\\ {\\bf 118} (1999) 213.\n";
  if ( !stratref.empty() ) ref() << stratref << '\n';
  for ( StringMap::iterator it = references.begin();
	it != references.end(); ++it )
    ref() << it->second << '\n';
  ref() << "\\end{thebibliography}\n"
	    << "\\end{document}" << endl;
  if ( !useStdout )
    reffile().close();
}

void EventGenerator::strategy(StrategyPtr s) {
  theStrategy = s;
}

int EventGenerator::count(const Exception & ex) {
  return ++theExceptions[make_pair(StringUtils::typeName(typeid(ex)),
				   ex.severity())];
}

void EventGenerator::printException(const Exception & ex) {
  switch ( ex.severity() ) {
  case Exception::info:
    log() << "* An information";
        break;
  case Exception::warning:
    log() << "* A warning";
    break;
  case Exception::setuperror:
    log() << "** A setup";
    break;
  case Exception::eventerror:
    log() << "** An event";
    break;
  case Exception::runerror:
    log() << "*** An run";
    break;
  case Exception::maybeabort:
  case Exception::abortnow:
    log() << "**** A serious";
    break;
  default:
    log() << "**** An unknown";
    break;
  }
  if ( ieve > 0 )
    log() << " exception of type " << StringUtils::typeName(typeid(ex))
	  << " occurred while generating event number "
	  << ieve << ": \n" << ex.message() << endl;
  else
    log() << " exception occurred in the initialization of "
	  << name() << ": \n" << ex.message() << endl;
  if ( ex.severity() == Exception::eventerror )
    log() << "The event will be discarded." << endl;
}

void EventGenerator::logWarning(const Exception & ex) {
  if ( ex.severity() != Exception::info &&
       ex.severity() != Exception::warning ) throw ex;
  ex.handle();  
  int c = count(ex);
  if ( c > maxWarnings ) return;
  printException(ex);
  if ( c == maxWarnings )
    log() << "No more warnings of this kind will be reported." << endl;
}
 
bool EventGenerator::
logException(const Exception & ex, tcEventPtr event) {
  bool noEvent = !event;
  ex.handle();
  int c = count(ex);
  if ( c <= maxWarnings ) {
    printException(ex);
    if ( c == maxWarnings )
      log() << "No more warnings of this kind will be reported." << endl;
  }
  if ( ex.severity() == Exception::info ||
       ex.severity() == Exception::warning ) {
    ex.handle();
    return false;
  }
  if ( ex.severity() == Exception::eventerror ) {
    if ( c < maxErrors || maxErrors <= 0 ) {
      ex.handle();
      if ( ThePEG_DEBUG_LEVEL > 0 && !noEvent ) log() << *event;
      return false;
    }
    if ( c > maxErrors ) printException(ex);
    log() << "Too many (" << c << ") exceptions of this kind has occurred. "
      "Execution will be stopped.\n";
  } else {
    log() << "This exception is too serious. Execution will be stopped.\n";
  }
  if ( !noEvent ) log() << *event;
  else log()
    << "An exception occurred before any event object was created!\n";
  dump();
  return true;
}

struct MatcherOrdering {
  bool operator()(tcPMPtr m1, tcPMPtr m2) {
    return m1->name() < m2->name() ||
      ( m1->name() == m2->name() && m1->fullName() < m2->fullName() );
  }
};

struct ObjectOrdering {
  bool operator()(tcIBPtr i1, tcIBPtr i2) {
    return i1->fullName() < i2->fullName();
  }
};

void EventGenerator::persistentOutput(PersistentOStream & os) const {
  set<tcPMPtr,MatcherOrdering> match(theMatchers.begin(), theMatchers.end());
  set<tcIBPtr,ObjectOrdering> usedset(usedObjects.begin(), usedObjects.end());
  os << theDefaultObjects << theLocalParticles << theStandardModel
     << theStrategy << theRandom << theEventHandler << theAnalysisHandlers
     << theHistogramFactory << theEventManipulator << thePath << theRunName
     << theNumberOfEvents << theObjectMap << theParticles
     << theQuickParticles << theQuickSize << match << usedset
     << ieve << weightSum << theDebugLevel << logNonDefault << printEvent
     << dumpPeriod << keepAllDumps << debugEvent
     << maxWarnings << maxErrors << theCurrentEventHandler
     << theCurrentStepHandler << useStdout << theIntermediateOutput << theMiscStream.str()
     << Repository::listReadDirs();
}

void EventGenerator::persistentInput(PersistentIStream & is, int) {
  string dummy;
  vector<string> readdirs;
  theGlobalLibraries = is.globalLibraries();
  is >> theDefaultObjects >> theLocalParticles >> theStandardModel
     >> theStrategy >> theRandom >> theEventHandler >> theAnalysisHandlers
     >> theHistogramFactory >> theEventManipulator >> thePath >> theRunName
     >> theNumberOfEvents >> theObjectMap >> theParticles
     >> theQuickParticles >> theQuickSize >> theMatchers >> usedObjects
     >> ieve >> weightSum >> theDebugLevel >> logNonDefault >> printEvent
     >> dumpPeriod >> keepAllDumps >> debugEvent
     >> maxWarnings >> maxErrors >> theCurrentEventHandler
     >> theCurrentStepHandler >> useStdout >> theIntermediateOutput >> dummy
     >> readdirs;
  theMiscStream.str(dummy);
  theMiscStream.seekp(0, std::ios::end);
  theObjects.clear();
  for ( ObjectMap::iterator it = theObjectMap.begin();
	it != theObjectMap.end(); ++it ) theObjects.insert(it->second);
  Repository::appendReadDir(readdirs);
}

void EventGenerator::setLocalParticles(PDPtr pd, int) {
  localParticles()[pd->id()] = pd;
}
  
void EventGenerator::insLocalParticles(PDPtr pd, int) {
  localParticles()[pd->id()] = pd;
}
  
void EventGenerator::delLocalParticles(int place) {
  ParticleMap::iterator it = localParticles().begin();
  while ( place-- && it != localParticles().end() ) ++it;
  if ( it != localParticles().end() ) localParticles().erase(it);
}

vector<PDPtr> EventGenerator::getLocalParticles() const {
  vector<PDPtr> ret;
  for ( ParticleMap::const_iterator it = localParticles().begin();
 	it != localParticles().end(); ++it ) ret.push_back(it->second);
  return ret;
}

void EventGenerator::setPath(string newPath) {
  if ( std::system(("mkdir -p " + newPath).c_str()) ) throw EGNoPath(newPath);
  if ( std::system(("touch " + newPath + "/.ThePEG").c_str()) )
    throw EGNoPath(newPath);
  if ( std::system(("rm -f " + newPath + "/.ThePEG").c_str()) )
    throw EGNoPath(newPath);
  thePath = newPath;
}

string EventGenerator::defPath() const {
  char * env = std::getenv("ThePEG_RUN_DIR");
  if ( env ) return string(env);
  return string(".");
}

ostream & EventGenerator::out() {
  return theOutStream;
}

ostream & EventGenerator::log() {
  return logfile().is_open()? logfile(): BaseRepository::cout();
}

ostream & EventGenerator::ref() {
  return reffile().is_open()? reffile(): BaseRepository::cout();
}

string EventGenerator::doSaveRun(string runname) {
  runname = StringUtils::car(runname);
  if ( runname.empty() ) runname = theRunName;
  if ( runname.empty() ) runname = name();
  EGPtr eg = Repository::makeRun(this, runname);
  string file =  eg->filename() + ".run";
  PersistentOStream os(file);
  os << eg;
  if ( !os ) return "Error: Save failed! (I/O error)";
  return "";
}

string EventGenerator::doMakeRun(string runname) {
  runname = StringUtils::car(runname);
  if ( runname.empty() ) runname = theRunName;
  if ( runname.empty() ) runname = name();
  Repository::makeRun(this, runname);
  return "";
}

bool EventGenerator::preinitRegister(IPtr obj, string fullname) {
  if ( !preinitializing ) throw InitException()
    << "Tried to register a new object in the initialization of an "
    << "EventGenerator outside of the pre-initialization face. "
    << "The preinitRegister() can only be called from a doinit() function "
    << "in an object for which preInitialize() returns true.";
  if ( objectMap().find(fullname) != objectMap().end() ) return false;
  obj->name(fullname);
  objectMap()[fullname] = obj;
  objects().insert(obj);
  obj->theGenerator = this;
  PDPtr pd = dynamic_ptr_cast<PDPtr>(obj);
  if ( pd ) theParticles[pd->id()] = pd;
  PMPtr pm = dynamic_ptr_cast<PMPtr>(obj);
  if ( pm ) theMatchers.insert(pm);
  return true;
}

IPtr EventGenerator::
preinitCreate(string classname, string fullname, string libraries) {
  if ( !preinitializing ) throw InitException()
    << "Tried to create a new object in the initialization of an "
    << "EventGenerator outside of the pre-initialization face. "
    << "The preinitCreate() can only be called from a doinit() function "
    << "in an object for which preInitialize() returns true.";
  if ( objectMap().find(fullname) != objectMap().end() ) return IPtr();
  const ClassDescriptionBase * db = DescriptionList::find(classname);
  while ( !db && libraries.length() ) {
    string library = StringUtils::car(libraries);
    libraries = StringUtils::cdr(libraries);
    DynamicLoader::load(library);
    db = DescriptionList::find(classname);
  }
  if ( !db ) return IPtr();
  IPtr obj = dynamic_ptr_cast<IPtr>(db->create());
  if ( !obj ) return IPtr();
  if ( !preinitRegister(obj, fullname) ) return IPtr();
  return obj;
}

string EventGenerator::
preinitInterface(IPtr obj, string ifcname, string cmd, string value) {
  if ( !preinitializing ) throw InitException()
    << "Tried to manipulate an external object in the initialization of an "
    << "EventGenerator outside of the pre-initialization face. "
    << "The preinitSet() can only be called from a doinit() function "
    << "in an object for which preInitialize() returns true.";
  if ( !obj ) return "Error: No object found.";
  const InterfaceBase * ifc = Repository::FindInterface(obj, ifcname);
  if ( !ifc ) return "Error: No such interface found.";
  try {
    return ifc->exec(*obj, cmd, value);
  }
  catch ( const InterfaceException & ex) {
    ex.handle();
    return "Error: " + ex.message();
  }
}

string EventGenerator::
preinitInterface(IPtr obj, string ifcname, int index,
		 string cmd, string value) {
  ostringstream os;
  os << index;
  return preinitInterface(obj, ifcname, cmd, os.str() + " " + value);
}

string EventGenerator::
preinitInterface(string fullname, string ifcname, string cmd, string value) {
  return preinitInterface(getObject<Interfaced>(fullname), ifcname, cmd, value);
}

string EventGenerator::
preinitInterface(string fullname, string ifcname, int index,
		 string cmd, string value) {
  return preinitInterface(getObject<Interfaced>(fullname), ifcname, index,
			  cmd, value);
}

tDMPtr EventGenerator::findDecayMode(string tag) const {
  for ( ObjectSet::const_iterator it = objects().begin();
	it != objects().end(); ++it ) {
    tDMPtr dm = dynamic_ptr_cast<tDMPtr>(*it);
    if ( dm && dm->tag() == tag ) return dm;
  }
  return tDMPtr();
}

tDMPtr EventGenerator::preinitCreateDecayMode(string tag) {
  return constructDecayMode(tag);
}

DMPtr EventGenerator::constructDecayMode(string & tag) {
  DMPtr rdm;
  DMPtr adm;
  int level = 0;
  string::size_type end = 0;
  while ( end < tag.size() && ( tag[end] != ']' || level ) ) {
    switch ( tag[end++] ) {
    case '[':
      ++level;
      break;
    case ']':
      --level;
      break;
    }
  }
  rdm = findDecayMode(tag.substr(0,end));
  if ( rdm ) return rdm;

  string::size_type next = tag.find("->");
  if ( next == string::npos ) return rdm;
  if ( tag.find(';') == string::npos ) return rdm;
  tPDPtr pd = getObject<ParticleData>(tag.substr(0,next));
  if ( !pd ) pd = findParticle(tag.substr(0,next));
  if ( !pd ) return rdm;
  rdm = ptr_new<DMPtr>();
  rdm->parent(pd);
  if ( pd->CC() ) {
    adm = ptr_new<DMPtr>();
    adm->parent(pd->CC());
    rdm->theAntiPartner = adm;
    adm->theAntiPartner = rdm;
  }
  bool error = false;
  tag = tag.substr(next+2);
  tPDPtr lastprod;
  bool dolink = false;
  do {
    switch ( tag[0] ) {
    case '[':
      {
	tag = tag.substr(1);
	tDMPtr cdm = constructDecayMode(tag);
	if ( cdm ) rdm->addCascadeProduct(cdm);
	else error = true;
      } break;
    case '=':
      dolink = true;
    case ',':
    case ']':
      tag = tag.substr(1);
      break;
    case '?':
      {
	next = min(tag.find(','), tag.find(';'));
	tPMPtr pm = findMatcher(tag.substr(1,next-1));
	if ( pm ) rdm->addProductMatcher(pm);
	else error = true;
	tag = tag.substr(next);
      } break;
    case '!':
      {
	next = min(tag.find(','), tag.find(';'));
	tPDPtr pd = findParticle(tag.substr(1,next-1));
	if ( pd ) rdm->addExcluded(pd);
	else error = true;
	tag = tag.substr(next);
      } break;
    case '*':
      {
	next = min(tag.find(','), tag.find(';'));
	tPMPtr pm = findMatcher(tag.substr(1,next-1));
	if ( pm ) rdm->setWildMatcher(pm);
	else error = true;
	tag = tag.substr(next);
      } break;
    default:
      {
	next = min(tag.find('='), min(tag.find(','), tag.find(';')));
	tPDPtr pdp = findParticle(tag.substr(0,next));
	if ( pdp ) rdm->addProduct(pdp);
	else error = true;
	tag = tag.substr(next);
	if ( dolink && lastprod ) {
	  rdm->addLink(lastprod, pdp);
	  dolink = false;
	}
	lastprod = pdp;
      } break;
    }
  } while ( tag[0] != ';' && tag.size() );
  if ( tag[0] != ';' || error ) {
    return DMPtr();
  }

  tag = tag.substr(1);
  
  DMPtr ndm = findDecayMode(rdm->tag());
  if ( ndm ) return ndm;
  pd->addDecayMode(rdm);
  if ( !preinitRegister(rdm, pd->fullName() + "/" + rdm->tag()) )
    return DMPtr();
  if ( adm ) {
    preinitRegister(adm, pd->CC()->fullName() + "/" + adm->tag());
    rdm->CC(adm);
    adm->CC(rdm);
  }

  return rdm;
}

tPDPtr EventGenerator::findParticle(string pdgname) const {
  for ( ParticleMap::const_iterator it = particles().begin();
	it != particles().end(); ++it )
    if ( it->second->PDGName() == pdgname ) return it->second;
  return tPDPtr();
}

tPMPtr EventGenerator::findMatcher(string name) const {
  for ( MatcherSet::const_iterator it = matchers().begin();
	it != matchers().end(); ++it )
    if ( (**it).name() == name ) return *it;
  return tPMPtr();
}

ClassDescription<EventGenerator> EventGenerator::initEventGenerator;

void EventGenerator::Init() {
  
  static ClassDocumentation<EventGenerator> documentation
    ("This is the main class used to administer an event generation run. "
     "The actual generation of each event is handled by the assigned "
     "<interface>EventHandler</interface> object. When the event generator"
     "is properly set up it can be initialized with the command "
     "<interface>MakeRun</interface> and/or saved to a file with the command "
     "<interface>SaveRun</interface>. If saved to a file, the event generator "
     "can be read into another program to produce events. The file can also "
     "be read into the <tt>runThePEG</tt> program where a number of events "
     "determined by the parameter <interface>NumberOfEvents</interface> is "
     "generated with each event analysed by the list of assigned "
     "<interface>AnalysisHandlers</interface>.");

  static Reference<EventGenerator,StandardModelBase> interfaceStandardModel
    ("StandardModelParameters",
     "The ThePEG::StandardModelBase object to be used to access standard "
     "model parameters in this run.",
     &EventGenerator::theStandardModel, false, false, true, false);

  static Reference<EventGenerator,EventHandler> interfaceEventHandler
    ("EventHandler",
     "The ThePEG::EventHandler object to be used to generate the "
     "individual events in this run.",
     &EventGenerator::theEventHandler, false, false, true, false);

  static RefVector<EventGenerator,AnalysisHandler> interfaceAnalysisHandlers
    ("AnalysisHandlers",
     "ThePEG::AnalysisHandler objects to be used to analyze the produced "
     "events in this run.",
     &EventGenerator::theAnalysisHandlers, 0, true, false, true, false);

  static Reference<EventGenerator,FactoryBase> interfaceHistogramFactory
    ("HistogramFactory",
     "An associated factory object for handling histograms to be used by "
     "<interface>AnalysisHandlers</interface>.",
     &EventGenerator::theHistogramFactory, true, false, true, true, true);

  static Reference<EventGenerator,EventManipulator> interfaceEventManip
    ("EventManipulator",
     "An ThePEG::EventManipulator called each time the generation of an "
     "event is stopped. The ThePEG::EventManipulator object is able to "
     "manipulate the generated event, as opposed to an "
     "ThePEG::AnalysisHandler which may only look at the event.",
     &EventGenerator::theEventManipulator, true, false, true, true);

  static RefVector<EventGenerator,ParticleData> interfaceLocalParticles
    ("LocalParticles",
     "Special versions of ThePEG::ParticleData objects to be used "
     "in this run. Note that to delete an object, its number in the list "
     "should be given, rather than its id number.",
     0, 0, false, false, true, false,
     &EventGenerator::setLocalParticles, &EventGenerator::insLocalParticles,
     &EventGenerator::delLocalParticles, &EventGenerator::getLocalParticles);

  static RefVector<EventGenerator,Interfaced> interfaceDefaultObjects
    ("DefaultObjects",
     "A vector of pointers to default objects. In a ThePEG::Reference or "
     "ThePEG::RefVector interface with the defaultIfNull() flag set, if a "
     "null pointer is encountered this vector is gone through until an "
     "acceptable object is found in which case the null pointer is replaced "
     "by a pointer to this object.",
     &EventGenerator::theDefaultObjects, 0, true, false, true, false, false);

  static Reference<EventGenerator,Strategy> interfaceStrategy
    ("Strategy",
     "An ThePEG::Strategy with additional ThePEG::ParticleData objects to "
     "be used in this run.",
     &EventGenerator::theStrategy, false, false, true, true);

  static Reference<EventGenerator,RandomGenerator> interfaceRandomGenerator
    ("RandomNumberGenerator",
     "An ThePEG::RandomGenerator object which should typically interaface to "
     "a CLHEP Random object. This will be the default random number generator "
     "for the run, but individual objects may use their own random generator "
     "if they wish.",
     &EventGenerator::theRandom, true, false, true, false);

  static Parameter<EventGenerator,string> interfacePath
    ("Path",
     "The directory where the output files are put.",
     &EventGenerator::thePath, ".", true, false,
      &EventGenerator::setPath, 0, &EventGenerator::defPath);
  interfacePath.directoryType();

  static Parameter<EventGenerator,string> interfaceRunName
    ("RunName",
     "The name of this run. This name will be used in the output filenames. "
     "The files wil be placed in the directory specified by the "
     "<interface>Path</interface> parameter"
     "If empty the name of the event generator will be used instead.",
     &EventGenerator::theRunName, "", true, false,
     0, 0, &EventGenerator::name);

  static Parameter<EventGenerator,long> interfaceNumberOfEvents
    ("NumberOfEvents",
     "The number of events to be generated in this run. If less than zero, "
     "the number of events is unlimited",
     &EventGenerator::theNumberOfEvents, 1000, -1, Constants::MaxInt,
     true, false, Interface::lowerlim);

  static Parameter<EventGenerator,int> interfaceDebugLevel
    ("DebugLevel",
     "The level of debug information sent out to the log file in the run. "
     "Level 0 only gives a limited ammount of warnings and error messages. "
     "Level 1 will print the first few events. "
     "Level 5 will print every event. "
     "Level 9 will print every step in every event.",
     &EventGenerator::theDebugLevel, 0, 0, 9, true, false, true);

  static Parameter<EventGenerator,int> interfacePrintEvent
    ("PrintEvent",
     "If the debug level is above zero, print the first 'PrintEvent' events.",
     &EventGenerator::printEvent, 0, 0, 1000, true, false, Interface::lowerlim);

  static Parameter<EventGenerator,long> interfaceDumpPeriod
    ("DumpPeriod",
     "If the debug level is above zero, dump the full state of the run every "
     "'DumpPeriod' events. Set it to -1 to disable dumping even in the case of errors.",
     &EventGenerator::dumpPeriod, 0, -1, Constants::MaxInt,
     true, false, Interface::lowerlim);

  static Switch<EventGenerator,bool> interfaceKeepAllDumps
    ("KeepAllDumps",
     "Whether all dump files should be kept, labelled by event number.",
     &EventGenerator::keepAllDumps, false, true, false);
  static SwitchOption interfaceKeepAllDumpsYes
    (interfaceKeepAllDumps,
     "Yes",
     "Keep all dump files, labelled by event number.",
     true);
  static SwitchOption interfaceKeepAllDumpsNo
    (interfaceKeepAllDumps,
     "No",
     "Keep only the latest dump file.",
     false);

  static Parameter<EventGenerator,long> interfaceDebugEvent
    ("DebugEvent",
     "If the debug level is above zero, step up to the highest debug level "
     "befor event number 'DebugEvent'.",
     &EventGenerator::debugEvent, 0, 0, Constants::MaxInt,
     true, false, Interface::lowerlim);

  static Parameter<EventGenerator,int> interfaceMaxWarnings
    ("MaxWarnings",
     "The maximum number of warnings of each type which will be printed.",
     &EventGenerator::maxWarnings,
     10, 1, 100, true, false, Interface::lowerlim);

  static Parameter<EventGenerator,int> interfaceMaxErrors
    ("MaxErrors",
     "The maximum number of errors of each type which will be tolerated. "
     "If more errors are reported, the run will be aborted.",
     &EventGenerator::maxErrors,
     10, -1, 100000, true, false, Interface::lowerlim);

  static Parameter<EventGenerator,long> interfaceQuickSize
    ("QuickSize",
     "The max absolute id number of particle data objects which are accessed "
     "quickly through a vector indexed by the id number.",
     &EventGenerator::theQuickSize,
     7000, 0, 50000, true, false, Interface::lowerlim);


  static Command<EventGenerator> interfaceSaveRun
    ("SaveRun",
     "Isolate, initialize and save this event generator to a file, from which "
     "it can be read in and run in another program. If an agument is given "
     "this is used as the run name, otherwise the run name is taken from the "
     "<interface>RunName</interface> parameter.",
     &EventGenerator::doSaveRun, true);


  static Command<EventGenerator> interfaceMakeRun
    ("MakeRun",
     "Isolate and initialize this event generator and give it a run name. "
     "If no argument is given, the run name is taken from the "
     "<interface>RunName</interface> parameter.",
     &EventGenerator::doMakeRun, true);

  interfaceEventHandler.rank(11.0);
  interfaceSaveRun.rank(10.0);
  interfaceMakeRun.rank(9.0);
  interfaceRunName.rank(8.0);
  interfaceNumberOfEvents.rank(7.0);
  interfaceAnalysisHandlers.rank(6.0);


  static Switch<EventGenerator,bool> interfaceUseStdout
    ("UseStdout",
     "Redirect the logging and output to stdout instead of files.",
     &EventGenerator::useStdout, false, true, false);
  static SwitchOption interfaceUseStdoutYes
    (interfaceUseStdout,
     "Yes",
     "Use stdout instead of log files.",
     true);
  static SwitchOption interfaceUseStdoutNo
    (interfaceUseStdout,
     "No",
     "Use log files.",
     false);
  

  static Switch<EventGenerator,int> interfaceLogNonDefault
    ("LogNonDefault",
     "Controls the printout of important interfaces which has been changed from their default values.",
     &EventGenerator::logNonDefault, -1, true, false);
  static SwitchOption interfaceLogNonDefaultYes
    (interfaceLogNonDefault,
     "Yes",
     "Always print changed interfaces.",
     1);
  static SwitchOption interfaceLogNonDefaultOnDebug
    (interfaceLogNonDefault,
     "OnDebug",
     "Only print changed interfaces if debugging is turned on.",
     0);
  static SwitchOption interfaceLogNonDefaultNo
    (interfaceLogNonDefault,
     "No",
     "Don't print changed interfaces.",
     -1);
  interfaceLogNonDefault.setHasDefault(false);

  static Switch<EventGenerator,bool> interfaceIntermediateOutput
    ("IntermediateOutput",
     "Modified event number count with the number of events processed so far, "
     "which updates at least every 10000 events, together with the corresponding "
     "intermediate estimate for the cross section plus the integration error.",
     &EventGenerator::theIntermediateOutput, false, true, false);
  static SwitchOption interfaceIntermediateOutputYes
    (interfaceIntermediateOutput,
     "Yes",
     "Show the modified event number count with the number of events processed so far, "
     "plus further information on the intermediate cross section estimate.",
     true);
  static SwitchOption interfaceIntermediateOutputNo
    (interfaceIntermediateOutput,
     "No",
     "Show the usual event number count with the number of events processed so far, "
     "but no further information on the intermediate cross section estimate.",
     false);

}

EGNoPath::EGNoPath(string path) {
  theMessage << "Cannot set the directory path for output files to '" << path
	     << "' because the directory did not exist and could not be "
	     << "created.";
  severity(warning);
}

