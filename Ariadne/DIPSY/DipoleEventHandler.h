// -*- C++ -*-
#ifndef DIPSY_DipoleEventHandler_H
#define DIPSY_DipoleEventHandler_H
//
// This is the declaration of the DipoleEventHandler class.
//

#include "ThePEG/Handlers/EventHandler.h"
#include "DipoleEventHandler.fh"
#include "WaveFunction.h"
#include "DipoleXSec.h"
#include "Emitter.h"
#include "Swinger.h"
#include "EventFiller.h"
// #include "DiffractiveEventFiller.h"
#include "ImpactParameterGenerator.h"
#include "DipoleAnalysisHandler.h"
#include "ThePEG/StandardModel/AlphaSBase.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Utilities/XSecStat.h"

namespace DIPSY { 

using namespace ThePEG;

/**
 * Here is the documentation of the DipoleEventHandler class.
 *
 * @see \ref DipoleEventHandlerInterfaces "The interfaces"
 * defined for DipoleEventHandler.
 */
class DipoleEventHandler: public EventHandler {

public:

  /** Declare a pointer to an AlphaSBase object. */
  typedef Ptr<AlphaSBase>::pointer ASPtr;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleEventHandler();

  /**
   * The destructor.
   */
  virtual ~DipoleEventHandler();
  //@}

public:

  /**
   * Initialize this event handler and all related objects needed to
   * generate events.
   */
  virtual void initialize();

  /**
   * Pre-sample the cross section before generating events.
   * This does inclusive cross section for non-diffractive collisions.
   * Diffractive inclusive cross sections are redirected.
   */
  void presample();

  /**
   * Pre-sample the single diffractive cross section before
   * generating events. This requires a different method, as it
   * needs a double loop over events.
   */
  void diffractivePresample();

  /**
   * Generate an event.
   */
  virtual EventPtr generateEvent();

  /**
   * Generate a diffractive event.
   */
  virtual EventPtr generateDiffractiveEvent();

  /**
   * Write out statistics.
   */
  virtual void statistics(ostream &) const;

  /**
   * The total integrated cross section of the processes generated in
   * this run. 
   * @return 0 if no integrated cross section could be estimated.
   */
  virtual CrossSection integratedXSec() const;

  /**
   * The error total integrated cross section of the processes
   * generated in this run.
   * @return 0 if no integrated cross section could be estimated.
   */
  virtual CrossSection integratedXSecErr() const;

  /**
   * The overestimated cross section ised in the generation.
   */
  CrossSection maxXSec() const;

  /**
   * Histogram scale. A histogram bin which has been filled with the
   * weights associated with the Event objects should be scaled by
   * this factor to give the correct cross section.
   */
  virtual CrossSection histogramScale() const;

  /**
   * Return the running coupling for the given size (scale).
   */
  double alphaS(InvEnergy r) const;

  /** @name Simple access functions. */
  //@{
  /**
   * The number of different colour indices.
   */
  inline int nColours() const {
    return theNColours;
  }

  /**
   * The general hadronic size.
   */
  inline InvEnergy rMax() const {
    return theRMax;
  }

  /**
   * Get the typical size of a baryon to be used by wave functions not
   * defining their own.
   */
  inline InvEnergy baryonSize() const {
    return theBaryonSize > ZERO? theBaryonSize: rMax();
  }

  /**
   * The maximum range at which coherent emissions are allowed.
   */
  inline InvEnergy coherenceRange() const {
    return theCoherenceRange;
  }

  /**
   * Return the effective parton mode.
   */
  inline int effectivePartonMode() const {
    return theEffectivePartonMode;
  }

  /**
   * Return the collision type.
   */
  inline int collisionType() const {
    return theCollisionType;
  }

  /**
   * The value of \f$\Lambda_{QCD}\f$ to be used in the running coupling.
   */
  inline Energy LambdaQCD() const {
    return theLambdaQCD;
  }

  /**
   * The number of flavours to be used in the running coupling.
   */
  inline int nF() const {
    return theNF;
  }

  /**
   * The value of the constant coupling. If zero, a running coupling is assumed.
   */
  inline double fixedAlphaS() const {
    return theFixedAlphaS;
  }

  /**
   * An external \f$\alpha_S\f$ object to be used if LambdaQCD() is zero.
   */
  inline ASPtr externalAlphaS() const {
    return theExternalAlphaS;
  }

  /**
   * Alpha bar = alphas*Nc/pi
   */
  inline double alphaBar(InvEnergy r) const {
    return alphaS(r)*3.0/M_PI;
  }

  /**
   * Get the wave function of the incoming particle along the positive z-axis.
   */
  inline WaveFunction & WFR() const {
    return *theWFR;
  }

  /**
   * Get the wave function of the incoming particle along the negative z-axis.
   */
  inline WaveFunction & WFL() const {
    return *theWFL;
  }

  /**
   * Get the object responsible for generating the impact parameters.
   */
  inline const ImpactParameterGenerator & bGen() const {
    return *theBGen;
  }

  /**
   * Indicate whether a fudge factor should be included to tame the
   * high-pt tail according to an approximate matrix element
   * correction.
   */
  inline bool fudgeME() const {
    return theFudgeME;
  }

  /**
   * Indicate in which frame the dipole systems should collide.
   */
  inline double yFrame() const {
    return theYFrame;
  }

  /**
   * Return the rapidity to use for the interaction frame.
   */
  double interactionFrame(double ymin, double ymax) const;  

  /**
   * Get the number of collisions to analyze in the presampling.
   */
  inline int preSamples() const {
    return thePreSamples;
  }

  /**
   * The number of left-moving systems to generate for each presample.
   */
  inline int preSampleL() const {
    return thePreSampleL;
  }

  /**
   * The number of right-moving systems to generate for each presample.
   */
  inline int preSampleR() const {
    return thePreSampleR;
  }

  /**
   * The number of impact parameters to generate for each presample.
   */
  inline int preSampleB() const {
    return thePreSampleB;
  }

  /**
   * Get the object responsible for calculating the cross section for
   * two colliding dipole systems.
   */
  inline const DipoleXSec & xSecFn() const {
    return *theXSecFn;
  }

  /**
   * Get the object responsible for generating and performing dipole
   * emissions of gluons.
   */
  inline const Emitter & emitter() const {
    return *theEmitter;
  }

  /**
   * Set the object responsible for generating and performing dipole
   * emissions of gluons.
   * Added by CF to access from emitter.
   */
  inline void emitter(EmitterPtr em) {
    theEmitter = em;
  }

  /**
   * Set the object responsible for generating and performing dipole
   * swings.
   * Added by CF to access from emitter.
   */
  inline void swinger(SwingerPtr sw) {
    theSwinger = sw;
  }

  /**
   * Get the object responsible for generating and performing dipole swings.
   */
  inline const Swinger & swinger() const {
    return *theSwinger;
  }

  /**
   * Get the object responsible for generating and performing dipole swings.
   */
  inline tSwingerPtr swingPtr() const {
    return theSwinger;
  }

  /**
   * Get the object responsible for filling an event with final state gluons.
   */
  inline const EventFiller & eventFiller() const {
    return *theEventFiller;
  }

  // /**
  //  * Get the object responsible for filling an event with final state gluons.
  //  */
  // inline const DiffractiveEventFiller & diffractiveEventFiller() const {
  //   return *theDiffractiveEventFiller;
  // };

  /**
   * adds an error to the counter
   **/
  inline void err() const {
    theNErr++;
  };

  /**
   * returns teh number of errors so far.
   **/
  inline const long nErr() const {
    return theNErr;
  };
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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The number of different colour indices.
   */
  int theNColours;

  /**
   * The general hadronic size.
   */
  InvEnergy theRMax;

  /**
   * The typical size of a baryon to be used by wave functions not
   * defining their own. If zero, theRMax will be used instead.
   */
  InvEnergy theBaryonSize;

  /**
   * The maximum range at which partons can emit coherently as effective partons.
   */
  InvEnergy theCoherenceRange;

  /**
   * The way the partons are grouped into effective partons.
   */
  int theEffectivePartonMode;

  /**
   * The type of collison. 0: non diffractive, 1: single diffractive
   */
  int theCollisionType;

  /**
   * if the history of every event should be studied manually.
   */
  bool doShowHistory;

  /**
   * The value of \f$\Lambda_{QCD}\f$ to be used in the running coupling.
   */
  Energy theLambdaQCD;

  /**
   * The number of flavours to be used in the running coupling.
   */
  int theNF;

  /**
   * The value of the constant coupling. If zero, a running coupling is assumed.
   */
  double theFixedAlphaS;

  /**
   * An external \f$\alpha_S\f$ object to be used if LambdaQCD() is zero.
   */
  ASPtr theExternalAlphaS;

  /**
   * The wave function of the incoming particle along the positive z-axis.
   */
  WaveFunctionPtr theWFR;

  /**
   * The wave function of the incoming particle along the negative z-axis.
   */
  WaveFunctionPtr theWFL;

  /**
   * The object responsible for generating the impact parameters.
   */
  ImpactParameterGeneratorPtr theBGen;

  /**
   * Indicate in which frame the dipole systems should collide.
   */
  double theYFrame;

  /**
   * Indicate whether a fudge factor should be included to tame the
   * high-pt tail according to an approximate matrix element correction.
   */
  bool theFudgeME;

  /**
   * The number of collisions to analyze in the presampling.
   */
  int thePreSamples;

  /**
   * The number of left-moving systems generated for each presample.
   */
  int thePreSampleL;

  /**
   * The number of right-moving systems generated for each presample.
   */
  int thePreSampleR;

  /**
   * The number of impact parameters generated for each presample.
   */
  int thePreSampleB;

  /**
   * The object responsible for calculating the cross section for two
   * colliding dipole systems.
   */
  DipoleXSecPtr theXSecFn;

  /**
   * The object responsible for generating and performing dipole
   * emissions of gluons.
   */
  EmitterPtr theEmitter;

  /**
   * The object responsible for generating and performing dipole swings.
   */
  SwingerPtr theSwinger;

  /**
   * The object responsible for filling an event with final state gluons.
   */
  EventFillerPtr theEventFiller;

  // /**
  //  * The object responsible for filling an event with diffractive final state gluons.
  //  */
  // DiffractiveEventFillerPtr theDiffractiveEventFiller;

  /**
   * A list of analysis objects to be applied in the presample phase.
   */
  vector<DipoleAnalysisHandlerPtr> analyses; 

  /**
   * Collect statistics.
   */
  XSecStat stats;

  /**
   * The number of error messages.
   **/
  mutable long theNErr;


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

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
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans) throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}

public:

  /** Return the elapsed number of seconds since lat call. */
  static double elapsed();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleEventHandler & operator=(const DipoleEventHandler &);

};

}

#endif /* DIPSY_DipoleEventHandler_H */
