// -*- C++ -*-
//
// StandardEventHandler.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_StandardEventHandler_H
#define ThePEG_StandardEventHandler_H
// This is the declaration of the StandardEventHandler class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Repository/Strategy.fh"
#include "ThePEG/Handlers/SamplerBase.fh"
#include "ThePEG/PDF/PartonBin.fh"
#include "ThePEG/MatrixElement/MEBase.fh"
#include "SubProcessHandler.fh"
#include "StandardXComb.fh"
#include "StandardEventHandler.fh"
#include "ThePEG/Utilities/XSecStat.h"
#include <fstream>

namespace ThePEG {

/**
 * The StandardEventHandler class is the main class for generating simple
 * events without overlayed collisions. It is derived from the
 * basic EventHandler class.
 *
 * Besides the standard doinit() method, the StandardEventHandler needs to be
 * separately initialized with the initialize() method. In the
 * dofinish() method statistics is written out to the EventGenerators
 * default output file.
 *
 * @see \ref StandardEventHandlerInterfaces "The interfaces"
 * defined for StandardEventHandler.
 * @see EventHandler
 * @see EventGenerator
 * @see Event
 * 
 */
class StandardEventHandler: public EventHandler {

public:

  /** A vector of <code>SubProcessHandler</code>s. */
  typedef vector<SubHdlPtr> SubHandlerList;

  /** A weighted list of pointers to StandardXComb objects. */
  typedef Selector<StdXCombPtr> XSelector;

  /** A vector of pointers to StandardXComb objects. */
  typedef vector<StdXCombPtr> XVector;

  /** A vector of cross sections. */
  typedef vector<CrossSection> XSVector;

  /** Map of pointers to StandardXComb objects indexed by pointers to
   *  the corresponding MEBase object. */
  typedef map<tMEPtr,XVector> MEXMap;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  StandardEventHandler();

  /**
   * Destructor.
   */
  virtual ~StandardEventHandler();
  //@}

public:

  /**
   * Initialize this event handler and all related objects needed to
   * generate events.
   */
  virtual void initialize();

  /**
   * Write out accumulated statistics about intergrated cross sections
   * and stuff.
   */
  virtual void statistics(ostream &) const;

  /**
   * Return the sampler assigned to this event handler.
   */
  tSamplerPtr sampler() { return theSampler; }

  /**
   * Return the sampler assigned to this event handler.
   */
  tcSamplerPtr sampler() const { return theSampler; }

  /**
   * Histogram scale. A histogram bin which has been filled with the
   * weights associated with the Event objects should be scaled by
   * this factor to give the correct cross section.
   */
  virtual CrossSection histogramScale() const;

  /**
   * The estimated total integrated cross section of the processes
   * generated in this run.
   * @return 0 if no integrated cross section could be estimated.
   */
  virtual CrossSection integratedXSec() const;

  /**
   * The estimated error int total integrated cross section of the
   * processes generated in this run. 
   * @return 0 if no integrated cross section error could be estimated.
   */
  virtual CrossSection integratedXSecErr() const;

  /**
   * The estimated total integrated cross section of the processes
   * generated in this run, excluding reweighting.
   * @return 0 if no integrated cross section could be estimated.
   */
  virtual CrossSection integratedXSecNoReweight() const;

  /**
   * The estimated error int total integrated cross section of the
   * processes generated in this run, excluding reweighting. 
   * @return 0 if no integrated cross section error could be estimated.
   */
  virtual CrossSection integratedXSecErrNoReweight() const;

  /** @name Functions used for the actual generation */
  //@{
  /**
   * Return the cross section for the chosen phase space point.
   * @param r a vector of random numbers to be used in the generation
   * of a phase space point.
   */
  virtual CrossSection dSigDR(const vector<double> & r);

  /**
   * Generate an event.
   */
  virtual EventPtr generateEvent();

  /**
   * Continue generating an event if the generation has been stopped
   * before finishing.
   */
  virtual EventPtr continueEvent();

  /**
   * Reweight a partially generated event.
   */
  void reweight(double factor) const;

  /**
   * Return the vector of StandardXComb objects.
   */
  const XVector & xCombs() const { return theXCombs; }

  /**
   * Change the XComb object
   */
  virtual void select(tXCombPtr newXComb);

  /**
   * Return the boost needed to transform the current event from the
   * CMS system to the lab system.
   */
  const LorentzRotation & currentEventBoost() const { return theCurrentEventBoost; }
  //@}

  /** @name Simple access functions */
  //@{
  /**
   * Return a reference to the Cuts of this
   * EventHandler. Note that these cuts may be overridden by the
   * SubProcess chosen.
   */
  tCutsPtr cuts() const { return theCuts; }

  /**
   * Return the number of separate bins of StandardXComb objects to
   * sample.
   */
  int nBins() const;

  /**
   * Return the number of phase space dimensions needed for the
   * sampling of indicated bin of StandardXComb objects.
   */
  int maxDim(int bin) const { return theMaxDims[bin]; }

  /**
   * The number of phase space dimensions used by the luminosity
   * function.
   */
  int lumiDim() const { return theLumiDim; }

  /**
   * The number of dimensions of the basic phase space to generate
   * sub-processes in for a given bin of StandardXComb objects.
   */
  int nDim(int bin) const { return lumiDim() + maxDim(bin); }
  //@}

protected:

  /**
   * Generate a phase space point and return the corresponding cross
   * section. Is called from sSigDR(const vector<double> &).
   * @param ll a pair of doubles giving the logarithms of the (inverse
   * energy fractions of the maximum CMS energy of the incoming
   * particles.
   * @param maxS the maximum squared CMS energy of the incoming particles.
   * @param ibin the preselected bin of StandardXComb objects to choose
   * sub-process from
   * @param nr the number of random numbers availiable in \a r.
   * @param r an array of random numbers to be used to generate a
   * phase-space point.
   */
  virtual CrossSection dSigDR(const pair<double,double> ll, Energy2 maxS,
		      int ibin, int nr, const double * r);

  /**
   * Select an StandardXComb. Given a preselected bin, \a ibin of
   * StandardXComb objects pick one to generate the corresponding
   * sub-process with the given \a weight.
   */
  tStdXCombPtr select(int bin, double & weight);

  /**
   * Create and add <code>StandardXComb</code> objects.
   *
   * @param maxEnergy the maximum CMS energy of the incoming particles.
   * @param sub a pointer to the SubProcessHandler object.
   * @param extractor a pointer to the PartonExtractor object.
   * @param cuts a pointer to the Cuts object.
   * @param ckkw a pointer to a CascadeHandler to be used for CKKW reweighting.
   * @param me a pointer to the MEBase object.
   * @param pBins a pair of <code>PartonBin</code>s describing the
   * partons extracted from the particles
   * @param allPBins all available parton bins at the given energy
   */
  void addME(Energy maxEnergy, tSubHdlPtr sub, tPExtrPtr extractor,
	     tCutsPtr cuts, tCascHdlPtr ckkw, tMEPtr me, const PBPair & pBins,
	     const PartonPairVec& allPBins);

  /**
   * For the sub-procss and phase-space point selected in the previous
   * call to dSigDR, produce the first step of an actual Collision.
   */
  tCollPtr performCollision();

  /**
   * Initialize groups of <code>StepHandler</code>s. This overrides
   * the method in the EventHandler, and the
   * <code>StepHandler</code>s given in the currently selected
   * SubProcess take precedence over the ones specified in the
   * EventHandler sub class.
   */
  virtual void initGroups();

  /**
   * Return the boost needed to transform the current collision from
   * the CMS system to the lab system. By default this is the unit
   * transformation, but an EventHandler derived from this class may
   * override it.
   */
  LorentzRotation & currentEventBoost() { return theCurrentEventBoost; }

  /**
   * Set information about the current sub-process.
   */
  void setScale(Energy2);

  /**
   * Return the vector of StandardXComb objects.
   */
  XVector & xCombs()  { return theXCombs; }

  /**
   * Throw away the last generated event before generating a new one.
   */
  virtual void clean();

private:

  /**
   * Access the list of sub-process handlers.
   */
  const SubHandlerList & subProcesses() const { return theSubProcesses; }

  /**
   * Access the list of sub-process handlers.
   */
  SubHandlerList & subProcesses() { return theSubProcesses; }

public:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  virtual void doupdate();

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
   * run has ended. Writes out statistics on the generation.
   */
  virtual void dofinish();
  //@}

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

  /**
   * Reject a (partially) generated event.
   * @param weight the weight given for the event.
   */
  void reject(double weight);

private:

  /**
   * The first of the incoming particle types.
   */

  PDPtr theIncomingA;
  /**
   * The second of the incoming particle types.
   */

  PDPtr theIncomingB;

  /**
   * The list of <code>SubProcessHandler</code>s.
   */
  SubHandlerList theSubProcesses;

  /**
   * The kinematical cuts used for this collision handler.
   */
  CutsPtr theCuts;

  /**
   * True if cuts on collision objects should be performed
   */
  bool collisionCuts;

  /**
   * The StandardXComb objects.
   */
  XVector theXCombs;

  /**
   * The number of degrees of freedom needed to generate the phase
   * space for the different bins.
   */
  vector<int> theMaxDims;

  /**
   * The boost needed to transform the current collision from the CMS
   * system to the lab system.
   */
  LorentzRotation theCurrentEventBoost;

  /**
   * The phase space sampler responsible for generating phase space
   * points according to the cross section given by this event
   * handler.
   */
  SamplerPtr theSampler;

  /**
   * The number of phase space dimensions used by the luminosity
   * function.
   */
  int theLumiDim;

  /**
   * The overall cross section statistics
   */
  mutable XSecStat xSecStats;

  /**
   * Standard Initialization object.
   */
  static ClassDescription<StandardEventHandler> initStandardEventHandler;

  /**
   * Helper function for the interface.
   */
  void setIncomingA(PDPtr);

  /**
   * Helper function for the interface.
   */
  void setIncomingB(PDPtr);

protected:

  /** @cond EXCEPTIONCLASSES */
  /**
   * Exception class used by EventHandler when a StepHandler of the
   * wrong class was added.
   */
  class StandardEventHandlerUpdateException: public UpdateException {};

  /**
   * Exception class used by EventHandler when a StepHandler of the
   * wrong class was added.
   */
  class StandardEventHandlerInitError: public Exception {};
  /** @endcond */

private:

  /**
   * Private and non-existent assignment operator.
   */
  const StandardEventHandler & operator=(const StandardEventHandler &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of StandardEventHandler.
 */
template <>
struct BaseClassTrait<StandardEventHandler,1>: public ClassTraitsType {
  /** Typedef of the base class of StandardEventHandler. */
  typedef EventHandler NthBase;
};

/**
 * The following template specialization informs ThePEG about the name
 * of theEventHandler class and the shared object where it is defined.
 */
template <>
struct ClassTraits<StandardEventHandler>
  : public ClassTraitsBase<StandardEventHandler> {
  /**
   * Return the class name.
   */
  static string className() { return "ThePEG::StandardEventHandler"; }
};

/** @endcond */

}

#endif /* ThePEG_StandardEventHandler_H */
