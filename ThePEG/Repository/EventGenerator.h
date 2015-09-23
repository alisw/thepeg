// -*- C++ -*-
//
// EventGenerator.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_EventGenerator_H
#define ThePEG_EventGenerator_H
// This is the declaration of the EventGenerator class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Utilities/Named.h"
#include "EventGenerator.fh"
#include "RandomGenerator.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/Strategy.h"
#include "ThePEG/Repository/CurrentGenerator.fh"
#include "ThePEG/Utilities/ClassDescription.h"
#include "ThePEG/Handlers/EventHandler.fh"
#include "ThePEG/Analysis/FactoryBase.fh"
#include <fstream>
#include "EventGenerator.xh"

namespace ThePEG {

/**
 * The EventGenerator class manages a whole event generator run. It
 * keeps a list of all Interfaced objects which are needed for a
 * particular run (these objects each have a pointer back to the
 * EventGenerator). Some objects are special, such as a default
 * RandomGenerator object, a StandardModelBase object and a Strategy
 * object and lists of ParticleData and MatcherBase objects used in
 * the run.
 *
 * The <code>EventGenerator</code> also manages information about the
 * run such as the exceptions being thrown, files to write output and
 * error messages to, etc.
 *
 * There are three main external member functions:<BR>
 * go() generates a specified number of events and exits.<BR>
 *
 * shoot() generates one Event and returns it.<BR>
 *
 * generateEvent() takes an initial Step or a partially generated
 * Event as argument and generates subsequent steps defined in the
 * generator.<BR>
 *
 * doShoot() is a virtual function called by shoot() and may be
 * overridden in sub-classes.<BR>
 *
 * doGenrateEvent() is a virtual function called by generateEvent() and
 * may to be overridden in sub-classes.
 *
 * @see \ref EventGeneratorInterfaces "The interfaces"
 * defined for EventGenerator.
 * @see Interfaced
 * @see RandomGenerator
 * @see StandardModelBase
 * @see Strategy
 * @see ParticleData
 * @see Event
 * @see Step
 * @see FullEventGenerator
 * 
 */
class EventGenerator: public Interfaced {

  /** The Repository is a friend. */
  friend class Repository;

public:

  /** A map of integers giving the number of times an exception of the
   *  key type has been thrown. */
  //typedef map<const type_info *, int> ExceptionMap;
  //typedef map<Exception, int, ExceptionComparison > ExceptionMap;
  typedef map<pair<string, Exception::Severity>, int> ExceptionMap;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  EventGenerator();

  /**
   * Copy-constructor.
   */
  EventGenerator(const EventGenerator &);

  /**
   * Destructor.
   */
  virtual ~EventGenerator();
  //@}

public:

  /** @name Access special objects in the run. */
  //@{
  /**
   * Return a pointer to the standard model parameters.
   */
  tSMPtr standardModel() const { return theStandardModel; }

  /**
   * Return a pointer to the strategy object containing a set of
   * non-default particles to use.
   */
  tStrategyPtr strategy() const { return theStrategy; }

  /**
   * Get the currently active EventHandler.
   */
  tEHPtr currentEventHandler() const { return theCurrentEventHandler; }

  /**
   * Set the currently active EventHandler.
   */
  void currentEventHandler(tEHPtr eh) { theCurrentEventHandler = eh; }

  /**
   * Get the currently active step handler.
   */
  tStepHdlPtr currentStepHandler() const { return theCurrentStepHandler; }

  /**
   * Set the currently active step handler.
   */
  void currentStepHandler(tStepHdlPtr sh) { theCurrentStepHandler = sh; }

  /**
   * Return a pointer to the EventHandler.
   */
  tEHPtr eventHandler() const { return theEventHandler; }

  /**
   * Return the vector of analysis objects to be used in the run.
   */
  AnalysisVector & analysisHandlers() { return theAnalysisHandlers; }

  /**
   * Return a pointer to an associated factory objects for handling
   * histograms to be used by <code>AnalysisHandler</code>s.
   */
  tHistFacPtr histogramFactory() const { return theHistogramFactory; }

  /**
   * Return the EventManipulator used in the run.
   */
  tEvtManipPtr manipulator() const { return theEventManipulator; }
  //@}

public:

  /** @name Main functions to controll the run. */
  //@{
  /**
   * Initialize this generator. This is done automatically if 'go()'
   * is used. Calls the virtual method doInitialize().
   */
  void initialize(bool initOnly = false);

  /**
   * Run this EventGenerator session. Calls the virtual method doGo().
   *
   * @param next the number of the firts event to be generated. If
   * negative it is assumed that this generator was previously
   * interrupted (or dumped to a file) and the execution will resume
   * from where it started. Default is 1.
   * @param maxevent the maximum number of events to be generated. If negative
   * the N() is used instead. Default is -1.
   * @param tics if true information the number of events generated
   * and elapsed time will be written to std::cerr after each event.
   */
  void go(long next = 1, long maxevent = -1, bool tics = false);

  /**
   * Generate one event. Calls the virtual method doShoot();
   */
  EventPtr shoot();

  /**
   * Finish generating an \a event which has already been partially
   * constructed from the outside.  Calls the virtual method do
   * doGenerateEvent().
   */
  EventPtr generateEvent(Event & event);

  /**
   * Finish generating an event starting from a \a step which has
   * already been partially constructed from the outside.  Calls the
   * virtual method do doGenerateEvent().
   */
  EventPtr generateEvent(Step & step);

  /**
   * Indicate that the run has ended and call finish() for all objects
   * including this one. Note that finish() should not be called
   * directly.
   */
  void finalize();

  /**
   * Dynamically load the Main class in the given \a file, making it
   * run its Init() method where it may use this EventGenerator. Also
   * call the initialize function before and the finish() function
   * afterwards.
   */
  bool loadMain(string file);

  /**
   * Return the maximum center of mass energy possible for an
   * event. Return zero if the assigned EventHander is not able to
   * generatr full events.
   */
  virtual Energy maximumCMEnergy() const;

  /**
   * The number of the event currently being generated.
   */
  long currentEventNumber() const { return ieve; }

  /**
   * Return the event being generated.
   */
  tcEventPtr currentEvent() const;

  /**
   * Dump the full state of the current run - including the number of
   * generated events, so that it can be fully continued from this point.
   */
  virtual void dump() const;

  /**
   * Register a given object as used. Only objects registered in this
   * way will be included in the file with model references.
   */
  void use(const Interfaced & i);

  /**
   * Set the random seed for the global random number generator. Also
   * set the interfaced member variable.
   */
  void setSeed(long seed);

  /**
   * Log a given exception.
   */
  void logWarning(const Exception &);

  /**
   * The number of events to be generated in this run.
   */
  long N() const { return theNumberOfEvents; }

  /**
   * Histogram scale. A histogram bin which has been filled with the
   * weights associated with the Event objects should be scaled by
   * this factor to give the correct cross section.
   */
  CrossSection histogramScale() const;

  /**
   * The total integrated cross section of the processes generated in
   * this run.
   */
  CrossSection integratedXSec() const;

  /**
   * The error estimate for the total integrated cross section of the
   * processes generated in this run.
   */
  CrossSection integratedXSecErr() const;

  /**
   * The sum of all weight of the events generated so far.
   */
  double sumWeights() const { return weightSum; }
  //@}

  /** @name Functions for accessing output files. */
  //@{
  /**
   * The base filename used in this run. The actual files are called
   * <code>filename.run</code>, <code>filename.dump</code>,
   * <code>filename.out</code>, <code>filename.log</code> and
   * <code>filename.tex</code> for the input configuration file,
   * output dump file, output file, log file, and reference
   * file respectively. The filename is constructed from the path()
   * and runName().
   */
  string filename() const { return path() + "/" + runName(); }

  /**
   * Return the name assigned to this run. If no name is given, the
   * name of the EventGenerator object is returned.
   */
  string runName() const { return theRunName.size()? theRunName: name(); }

  /**
   * The directory in which the filename() is located
   */
  string path() const { return thePath; }

  /**
   * Has the generator been asked to redirect everything to standard
   * output?
   */
  bool useStdOut() const { return useStdout; }

  /**
   * Open all ouput files.
   */
  void openOutputFiles();

  /**
   * Flush the content of the internal output string stream to the .out file.
   */
  void flushOutputFile();

  /**
   * Close all ouput files.
   */
  void closeOutputFiles();

  /**
   * Return a reference to the output file stream.
   */
  ofstream & outfile() { return theOutfile; }

  /**
   * Return a reference to the log file stream.
   */
  ofstream & logfile() { return theLogfile; }

  /**
   * Return a reference to the reference file stream. This file is
   * used to output LaTeX text with information about the models used
   * in the run.
   */
  ofstream & reffile() { return theReffile; }

  /**
   * This stream should be used for output of information and
   * statistics of an EventGenerator run in the finish() phase, after
   * the actual generation has finished. When used at other times, the
   * output will be cashed internally before written out in the
   * finish() phase. This is then written to the .out file, or if
   * useStdOut() is true, to BaseRepository::cout().
   */
  ostream & out();

  /**
   * Return a reference to the stream connected to the file for logging
   * information. If no file is connected, BaseRepository::cout() will
   * be used instead.
   */
  ostream & log();

  /**
   * Return a reference to a stream to be used to redirect cout for
   * external modules which prints out messages there. The output will
   * instead be appended to the log() stream at the end of the run.
   */
  ostream & misc() {
    return theMiscStream;
  }

  /**
   * Return a reference to the stream connected to the filea for
   * references from used objects. If no file is connected,
   * BaseRepository::cout() will be used instead.
   */
  ostream & ref();
  //@}

  /** @name Access objects included in this run. */
  //@{
  /**
   * Return the set of objects used in this run.
   */
  const ObjectSet & objects() const { return theObjects; }


  /**
   * Return the map of objects used in this run indexed by their name.
   */
  const ObjectMap & objectMap() const { return theObjectMap; }

  /**
   * Return a garbage collected pointer to a given object. If the
   * object is not included in the run, a null pointer will be
   * returned.
   */
  template <typename T>
  typename Ptr<T>::pointer getPtr(const T &) const;

  /**
   * Return a pointer to an object present in this run given its full
   * name. Return the null pointer if non-existent.
   */
  IBPtr getPointer(string name) const;

  /**
   * Return a pointer to an object of type T present in this run given
   * its full name. Return the null pointer if non-existent. Calls
   * getPointer(string) and dynamically casts the result to the
   * requested pointer type.
   */
  template <typename T>
  typename Ptr<T>::pointer getObject(string name) const {
    return dynamic_ptr_cast<typename Ptr<T>::pointer>(getPointer(name));
  }

  /**
   * Return the default object for class T. Returns the null pointer
   * if non-existent.
   */
  template <typename T>
  typename Ptr<T>::pointer getDefault() const;

  /**
   * Create a particle instance corresponding to the given \a id
   * number.
   */
  PPtr getParticle(PID id) const;

  /**
   * Return a pointer to the ParticleData object corresponding to the
   * given \a id number.
   */
  PDPtr getParticleData(PID id) const;

  /**
   * Return a reference to the complete list of matchers in this
   * generator.
   */
  const MatcherSet & matchers() const { return theMatchers; }

  /**
   * Return a reference to the complete map of particle data objects
   * in this generator, indexed by their id numbers.
   */
  const ParticleMap & particles() const { return theParticles; }

  /**
   * Return a reference to the set of objects which have been
   * registered as used during the current run.
   */
  const ObjectSet & used() const { return usedObjects; }
  //@}

protected:

  /**
   * Check if there has been an interrupt signal from the OS.
   * If that's the case, finalize() is called
   */
  void checkSignalState();

  /**
   * Return a reference to the default RandomGenerator object in this
   * run.
   */
  RandomGenerator & random() const { return *theRandom; }

  /**
   * Finish the setup of an event generator run. Set run name, all
   * particles, matchers and other objects to be used. Is used by the
   * Repository when isolating an EventGenerator.
   */
  void setup(string newRunName, ObjectSet & newObjects,
	     ParticleMap & newParticles, MatcherSet & newMatchers);

  /** @name Main virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Run this EventGenerator session. Is called from go(long,long,bool).
   */
  virtual void doGo(long next, long maxevent, bool tics);

  /**
   * Initialize this generator. Is called from initialize().
   */
  virtual void doInitialize(bool initOnly = false);

  /**
   * Generate one event. Is called from shoot().
   */
  virtual EventPtr doShoot();

  /**
   * Write out the number of events generated and the elapsed time in
   * suitable periods.
   */
  void tic(long currev = 0, long totev = 0) const;

  /**
   * Finish generating an event constructed from the outside. Is
   * called by generateEvent(tEventPtr).
   */
  virtual EventPtr doGenerateEvent(tEventPtr);

  /**
   * Finish generating an event starting from a Step constructed from
   * the outside. Is called by generateEvent(tStepPtr).
   */
  virtual EventPtr doGenerateEvent(tStepPtr);
  //@}

  /**
   * Print the message of an exception to the log file.
   */
  void printException(const Exception &);

  /**
   * Log a given exception.
   */
  bool logException(const Exception &, tcEventPtr);

  /**
   * Set number of events to be generated.
   */
  void N(long n) { theNumberOfEvents = n; }

  /**
   * Set the name of this run
   */
  void runName(string f) { theRunName = f; }

public:

  /**
   * Append a tag to the run name. Derived classes may put special
   * meaning to the tags. 
   */
  virtual void addTag(string tag) {
    runName(runName() + tag);
  }

private:

  /**
   * Return the vector of default objects.
   */
  const vector<IPtr> & defaultObjects() const { return theDefaultObjects; }

  /**
   * Access the special particles used in this generator. Not relevant
   * in the run phase.
   */
  ParticleMap & localParticles() { return theLocalParticles; }

  /**
   * Access the special particles used in this generator. Not relevant
   * in the run phase.
   */
  const ParticleMap & localParticles() const { return theLocalParticles; }

  /**
   * Set the directory where the output files will be stored.
   */
  void path(string f) { thePath = f; }

  /**
   * Set a pointer to the strategy object containing a set of
   * non-default particles to use.
   */
  void strategy(StrategyPtr);

  /**
   * Isolate, initialize and save this generator to a file.
   */
  string doSaveRun(string);

  /**
   * Isolate and initialize this generator.
   */
  string doMakeRun(string);

public:

  /** @name The following functions may be called by objects belonging
      to this event generator during the initialization phase (in the
      doinit() function). It is typically used by objects which need
      to introduce other Interfaced objects depending the parameters
      of the StandardModel object used. Note that objects which use
      these functions <b>MUST</b> override the preInitialize()
      function to return true, otherwize the whole initialization
      procedure may be corrupted. */
  //@{
  /**
   * Register a new object to be included in the run currently being
   * initialized.
   *
   * @param obj (pointer to) the object being registered.
   *
   * @param fullname the full name including the directory path. Note
   * that although the full path is given the object will not be
   * inserted in the Repository, but only in this current
   * EventGenerator.
   *
   * @return false if another object of that name already exists.
   */
  bool preinitRegister(IPtr obj, string fullname);

  /**
   * Create a new Interfaced object to be used in the run being
   * initialized.
   *
   * @param classname the class name of the object being created.
   *
   * @param fullname the full name including the directory path. Note
   * that although the full path is given the object will not be
   * inserted in the Repository, but only in this current
   * EventGenerator.
   *
   * @param libraries an optional list of shared libraries to be
   * loaded to be able to create an object of the specified class.
   *
   * @return the created object if the it was successfully
   * created. Return null if the object could not be created or if
   * another object of that name already exists.
   */
  IPtr preinitCreate(string classname, string fullname,	string libraries = "");


  /**
   * Manipulate an interface of an Interfaced object.
   *
   * @param fullname the name including the full path of an object to
   * be manipulated.
   *
   * @param ifcname the name of the interface to be used.
   *
   * @param cmd the operation to be performed on the interface (set or
   * get).
   *
   * @param value Optional value to be passed to the interface.
   *
   * @return a string containing the result of the operation. If this
   * string starts with "Error: " then something went wrong.
   */
  string preinitInterface(string fullname, string ifcname, string cmd,
			  string value);

  /**
   * Manipulate an interface of vector type (RefVector or ParVector)
   * of an Interfaced object.
   *
   * @param fullname the name including the full path of an object to
   * be manipulated.
   *
   * @param ifcname the name of the interface to be used.
   *
   * @param index the vector index corresponding to the element to be
   * manipulated.
   *
   * @param cmd the operation to be performed on the interface (set,
   * get, insert or erase).
   *
   * @param value Optional value to be passed to the interface.
   *
   * @return a string containing the result of the operation. If this
   * string starts with "Error: " then something went wrong.
   */
  string preinitInterface(string fullname, string ifcname, int index,
			  string cmd, string value);

  /**
   * Manipulate an interface of an Interfaced object.
   *
   * @param obj the object to be manipulated.
   *
   * @param ifcname the name of the interface to be used.
   *
   * @param cmd the operation to be performed on the interface (set or
   * get).
   *
   * @param value Optional value to be passed to the interface.
   *
   * @return a string containing the result of the operation. If this
   * string starts with "Error: " then something went wrong.
   */
  string preinitInterface(IPtr obj, string ifcname, string cmd, string value);

  /**
   * Manipulate an interface of vector type (RefVector or ParVector)
   * of an Interfaced object.
   *
   * @param obj the object to be manipulated.
   *
   * @param ifcname the name of the interface to be used.
   *
   * @param index the vector index corresponding to the element to be
   * manipulated.
   *
   * @param cmd the operation to be performed on the interface (set,
   * get, insert or erase).
   *
   * @param value Optional value to be passed to the interface.
   *
   * @return a string containing the result of the operation. If this
   * string starts with "Error: " then something went wrong.
   */
  string preinitInterface(IPtr obj, string ifcname, int index,
			  string cmd, string value);

  /**
   * Find a decaymode given a decay \a tag.
   * @return null if no decay mode was found.
   */
  tDMPtr findDecayMode(string tag) const;

  /**
   * Create a decay mode according to the given tag.
   * @return null if no decay mode could be created.
   */
  tDMPtr preinitCreateDecayMode(string tag);

  /**
   * Find a particle in this run, using its PDG name.
   * @return null if no particle is found.
   */
  tPDPtr findParticle(string pdgname) const;

  /**
   * Find a matcher in this run given its \a name.
   * @return null if no mather is found.
   */
  tPMPtr findMatcher(string name) const;

private:

  /**
   * Used internally by preinitCreateDecayMode();
   */
  DMPtr constructDecayMode(string & tag);

  //@}

public:


  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);

  /**
   * The global libraries needed for objects used in this EventGenerator.
   */
  const vector<string> & globalLibraries() const {
    return theGlobalLibraries;
  }

  //@}

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

  /**
   * Additional things to do at the very end after the (do)finish(),
   * such as closing output files etc.
   */
  void finally();

  //@}

  /**
   * Return the set of all objects to be used in this run.
   */
  ObjectSet & objects() { return theObjects; }

  /**
   * Return the map of all objects to be used in this run indexed by
   * their name.
   */
  ObjectMap & objectMap() { return theObjectMap; }

  /**
   * Print out the .tex file with descriptions of and references to
   * all models used in the run.
   */
  void generateReferences();

  /**
   * Increase and return the count for the given exception.
   */
  int count(const Exception &);

private:


  /**
   * A vector of default objects.
   */
  vector<IPtr> theDefaultObjects;

  /**
   * Map of non-default particles used in this EventGenerator.
   */
  ParticleMap theLocalParticles;

  /**
   * Pointer to an object containing standard model parameters.
   */
  SMPtr theStandardModel;

  /**
   * Pointer to a strategy object with other non-default particles to
   * be used in this EventGenerator.
   */
  StrategyPtr theStrategy;

  /**
   * Pointer to the default RandomGenerator to be used in this run.
   */
  RanGenPtr theRandom;

  /**
   * Pointer to the event handler used to generate the indivudual
   * events.
   */
  EHPtr theEventHandler;

  /**
   * A vector of all analysis handlers to be called after each event.
   */
  AnalysisVector theAnalysisHandlers;

  /**
   * A pointer to an associated factory objects for handling
   * histograms to be used by <code>AnalysisHandler</code>s.
   */
  HistFacPtr theHistogramFactory;

  /**
   * A pointer to an optional event manipulator object.
   */
  EvtManipPtr theEventManipulator;

  /**
   * The directory where the input and output files resides.
   */
  string thePath;

  /**
   * The name of this run.
   */
  string theRunName;

  /**
   * A reference to the output file stream.
   */
  ofstream theOutfile;

  /**
   * A reference to the log file stream.
   */
  ofstream theLogfile;

  /**
   * A reference to the reference file stream.
   */
  ofstream theReffile;

  /**
   * A stream to be used to redirect cout for external modules which
   * prints out messages there. The output will instead be appended to
   * the log() stream at the end of the run.
   */
  ostringstream theMiscStream;

  /**
   * A string stream used as a buffer for messages written to the .out
   * file. The .out file should in rinciple only be written to in the
   * end of a run, during the finish() phase, but if anything is
   * written before that, it will be cashed in this string stream
   * before written out properly in the end of the run.
   */
  ostringstream theOutStream;

  /**
   * Remember the name of the file where the output should be
   * sent. This is set int openOutputFiles().
   */
  string theOutFileName;

  /**
   * Number of events to be generated in this run.
   */
  long theNumberOfEvents;

  /**
   * The set of all objects to be used in this run.
   */
  ObjectSet theObjects;

  /**
   * All objects to be used in this run mapped to their name.
   */
  ObjectMap theObjectMap;

  /**
   * The map of all particles to be used in this run, indexed by the
   * id number.
   */
  ParticleMap theParticles;
  /**
   * A vector of particles indexed by the id number for quick access.
   * Only particles with id number less than theQuickSize are
   * available.
   */
  PDVector theQuickParticles;

  /**
   * Only particles with id number less than theQuickSize are
   * available in theQuickParticles.
   */
  long theQuickSize;

  /**
   * A flag to tell if we are in the pre-initialization phase where
   * objects with preInitialize() functions returning true are
   * initialized before others.
   */
  bool preinitializing;

  /**
   * The set of all matchers to be used in this run.
   */
  MatcherSet theMatchers;

  /**
   * The set of objects which have actually been used in this run.
   */
  ObjectSet usedObjects;

protected:

  /**
   * The current event number;
   */
  long ieve;

  /**
   * The sum of the weights of the events produced so far.
   */
  double weightSum;

  /**
   * The debug level.
   */
  int theDebugLevel;

private:

  /**
   * List all modified interfaces in the log file. If positive always
   * do this, if negative never do it. If zero, only do it if
   * debugging is turned on.
   */
  int logNonDefault;

  /**
   * If the debug level is higher than 0, print the first 'printEvent'
   * events to the logfile.
   */
  int printEvent;

  /**
   * If the debug level is higher than 0, dump the complete state of
   * this run to the default dump file every 'dumpPeriod' events.
   * If 'dumpPeriod' is -1, dumping is disabled completely,
   * even when runs are aborted.
   */
  long dumpPeriod;

  /**
   * If this flag is true, keep all dump files of the run, 
   * labelled by event number.
   */
  bool keepAllDumps;

  /**
   * If the debug level is higher than 0, step up to the highest debug
   * level just before the event with number debugEvent is performed.
   */
  long debugEvent;

  /**
   * The maximum number of warnings reported of each type. If more
   * than maxWarnings warnings of one type is issued, the generation
   * will continue without reporting this warning.
   */
  int maxWarnings;

  /**
   * The maximum number of warnings and errors reported of each
   * type. If more than maxErrors errors is reported for one type the
   * run will be aborted. Disable the check by setting to -1.
   */
  int maxErrors;

  /**
   * A map of all Exceptions which have been caught by the event
   * generator and the number of time each exception type has been
   * caught.
   */
  ExceptionMap theExceptions;

private:

  /**
   * Utility function for the interface.
   */
  void setLocalParticles(PDPtr pd, int);

  /**
   * Utility function for the interface.
   */
  void insLocalParticles(PDPtr pd, int);

  /**
   * Utility function for the interface.
   */
  void delLocalParticles(int place);

  /**
   * Utility function for the interface.
   */
  vector<PDPtr> getLocalParticles() const;

  /**
   * Utility function for the interface.
   */
  void setPath(string newPath);

  /**
   * Utility function for the interface.
   */
  string defPath() const;

  /**
   * The UseRandom object constructed for the duration of an
   * EventGenerator run so that the default random number generator
   * always can be accessed through the static methods of the
   * UseRandom class.
   */
  UseRandom * theCurrentRandom;

  /**
   * The CurrentGenerator object constructed for the duration of an
   * EventGenerator run so that the default event generator always can
   * be accessed through the static methods of the CurrentGenerator
   * class.
   */
  CurrentGenerator * theCurrentGenerator;

  /**
   * The currently active EventHandler.
   */
  tEHPtr theCurrentEventHandler;

  /**
   * The currently active step handler.
   */
  tStepHdlPtr theCurrentStepHandler;


  /**
   * Whether to use files or stdout for logging and output.
   */
  bool useStdout;

  /**
   * Whether to use a modified event number count.
   */
  bool theIntermediateOutput;

  /**
   * The global libraries needed for objects used in this EventGenerator.
   */
  vector<string> theGlobalLibraries;

private:

  /**
   * Describe an abstract class with persistent data.
   */
  static ClassDescription<EventGenerator> initEventGenerator;

  /**
   *  Private and non-existent assignment operator.
   */
  EventGenerator & operator=(const EventGenerator &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of EventGenerator. */
template <>
struct BaseClassTrait<EventGenerator,1>: public ClassTraitsType {
  /** Typedef of the first base class of EventGenerator. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  EventGenerator class. */
template <>
struct ClassTraits<EventGenerator>: public ClassTraitsBase<EventGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::EventGenerator"; }
};

/** @endcond */

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "EventGenerator.tcc"
#endif

#endif /* ThePEG_EventGenerator_H */
