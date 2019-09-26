// -*- C++ -*-
#ifndef THEPEG_ProgressLog_H
#define THEPEG_ProgressLog_H
//
// This is the declaration of the ProgressLog class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"

namespace ThePEG {

/**
 * The ProgressLog class will not perform an actual analysis. Instead
 * it will write out a progress status on the standard log file. By
 * default it will write on event 1, 2, 5, 10, 20, 50, ... etc. But
 * optionally it can in addition also write out every given number of
 * seconds.
 *
 * The status line which is written out contains the current date and
 * time, the number of events processed so far and the total number of
 * events to be generated, two estimates of the time of completion
 * (one based on the current cpu usage and one based on the average
 * cpu usage [the usage is given in brackets]), and the host on which
 * the program is running, together with its process number.
 *
 * @see \ref ProgressLogInterfaces "The interfaces"
 * defined for ProgressLog.
 */
class ProgressLog: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ProgressLog();

  /**
   * The destructor.
   */
  virtual ~ProgressLog();
  //@}

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);
  //@}

  /**
   * Return the cpu clock in seconds.
   */
  static double fclock();

  /**
   * Check if it is time to write out a status line.
   */
  bool statusTime(long i, long n) const;

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
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
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
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * If larger than 0, a status line will be written every secstep second.
   */
  int secstep;

  /**
   * The clock when the run was started.
   */
  time_t time0;

  /**
   * The cpu clock when the run was started.
   */
  double fcpu0;

  /**
   * The clock the last time a status line was written out.
   */
  time_t time1;

  /**
   * The cpu clock the last time a status line was written out.
   */
  double fcpu1;

  /**
   * The host on which we are running.
   */
  string host;

  /**
   * The pid of the current process.
   */
  pid_t pid;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ProgressLog> initProgressLog;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ProgressLog & operator=(const ProgressLog &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ProgressLog. */
template <>
struct BaseClassTrait<ProgressLog,1> {
  /** Typedef of the first base class of ProgressLog. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ProgressLog class and the shared object where it is defined. */
template <>
struct ClassTraits<ProgressLog>
  : public ClassTraitsBase<ProgressLog> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::ProgressLog"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ProgressLog is implemented. It may also include several, space-separated,
   * libraries if the class ProgressLog depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "ProgressLog.so"; }
};

/** @endcond */

}

#endif /* THEPEG_ProgressLog_H */
