// -*- C++ -*-
#ifndef DIPSY_PTAnalysis_H
#define DIPSY_PTAnalysis_H
//
// This is the declaration of the PTAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include <iostream>
#include <fstream>

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the PTAnalysis class.
 *
 * @see \ref PTAnalysisInterfaces "The interfaces"
 * defined for PTAnalysis.
 */
class PTAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  PTAnalysis();

  /**
   * The destructor.
   */
  virtual ~PTAnalysis();
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

  /**
   * Transform the event to the desired Lorentz frame and return the
   * corresponding LorentzRotation.
   * @param event a pointer to the Event to be transformed.
   * @return the LorentzRotation used in the transformation.
   */
  virtual LorentzRotation transform(tEventPtr event) const;

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  virtual void analyze(const tPVector & particles, double weight);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   */
  virtual void analyze(tPPtr particle, double weight);
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

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /** The pt-distribution of charged particles -7 < y < -5. */
  tH1DPtr hpt75;

  /** The pt-distribution of charged particles -5 < y < -3. */
  tH1DPtr hpt53;

  /** The pt-distribution of charged particles -3 < y < -1. */
  tH1DPtr hpt31;

  /** The pt-distribution of charged particles -1 < y <  1. */
  tH1DPtr hpt11, hptm1, hptg1, hptG1;

  /** The pt-distribution of charged particles  1 < y <  3. */
  tH1DPtr hpt13;

  /** The pt-distribution of charged particles  3 < y <  5. */
  tH1DPtr hpt35;

  /** The pt-distribution of charged particles  5 < y < 07. */
  tH1DPtr hpt57;

  /** The eta-distribution of charged particles pt > 5 GeV. */
  tH1DPtr heta5;

  /**
   * The sum of weights.
   */
  double sumw;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PTAnalysis & operator=(const PTAnalysis &);

};

}

#endif /* DIPSY_PTAnalysis_H */
