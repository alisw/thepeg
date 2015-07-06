// -*- C++ -*-
#ifndef DIPSY_AnalysisProgress_H
#define DIPSY_AnalysisProgress_H
//
// This is the declaration of the AnalysisProgress class.
//

#include "Ariadne/DIPSY/DipoleAnalysisHandler.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * AnalysisProgress writes in the log file the progress of the
 * presampling phase of the DipoleEventHandler.
 *
 * @see \ref AnalysisProgressInterfaces "The interfaces"
 * defined for AnalysisProgress.
 */
class AnalysisProgress: public DIPSY::DipoleAnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  AnalysisProgress();

  /**
   * The destructor.
   */
  virtual ~AnalysisProgress();
  //@}

public:

  /** @name Standard virtual functions inherited from the classes. */
  //@{
  /**
   * Initialize the analysis object.
   */
  virtual void initialize();

  /**
   * Just write out progress information
   */
  virtual void analyze(const vector<DipoleStatePtr> &, const vector<DipoleStatePtr> &,
	const vector<ImpactParameters> &, const DipoleXSec &,
	const Vec3D &, double);

  /**
   * Finalize the analysis, (compute statistics etc.). \a neve is the
   * number of times analyze() has been called since last
   * initialize().
   */
  virtual void finalize(long neve) {}
  //@}


  /**
   * Return the cpu clock in seconds.
   */
  static double fclock();

  /**
   * Check if it is time to write out a status line.
   */
  bool statusTime(long i, long n) const;

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

private:

  /**
   * The current event.
   */
  int ieve;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AnalysisProgress & operator=(const AnalysisProgress &);

};

}

#endif /* DIPSY_AnalysisProgress_H */
