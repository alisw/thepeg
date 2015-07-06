// -*- C++ -*-
#ifndef Ariadne5_AriadneHandler_H
#define Ariadne5_AriadneHandler_H
//
// This is the declaration of the AriadneHandler class.
//

#include "Ariadne/Config/Ariadne5.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/StandardModel/AlphaSBase.h"
#include "AriadneHandler.fh"
#include "History.fh"
#include "ReweightBase.fh"
#include "DipoleState.fh"
#include "EmitterBase.fh"
#include "BornCheckerBase.fh"
#include "ScaleSetter.fh"
#include "DISFinder.fh"
#include "QCDDipoleFinder.fh"
#include "EMDipoleFinder.fh"
#include "ResonanceFinder.fh"
#include "History.h"
#include "ConsistencyChecker.fh"
#include "QCDDipole.fh"
#include "EMDipole.fh"
#include "Ariadne/Cascade/Models/RemnantModel.fh"
#include "Ariadne/Cascade/Models/ColourResonanceModel.fh"
#include "ThePEG/Analysis/FactoryBase.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The AriadneHandler class inherits form the ThePEG::CascadeHandler
 * base class and implements the Dipole Cascade Model for partonic
 * cascades.
 *
 * @see \ref AriadneHandlerInterfaces "The interfaces"
 * defined for AriadneHandler.
 */
class AriadneHandler: public CascadeHandler {

public:

  /**
   * Enum different options for running \f$\alpha_s\f$.
   */
  enum RunningOption {
    externalRunning = -1, /**< Use the \f$\alpha_s\f$ specified in the
			       current StandardModel object. */
    noRunning = 0,        /**< Use a fixed \f$\alpha_s\f$ specified by
			       alpha0(). */
    simpleRunning = 1,    /**< Use simple leading order running with
			       \f$\Lambda_{QCD}\f$ given by lambdaQCD() */
    internalRunning = 2   /**< Use internal AlphaSBase object in
			       theInternalAlphaS for the running. */
  };

  /**
   * Enumerate strategies for purging gluons with too small transverse momentum.
   */
  enum PurgeStrategy {
    onlybefore = -1,     /**<   Only purge gluons before the cascade. */
    onlyafter = -2,      /**<   Only purge gluons after the cascade. */
    neverpurge = 0,           /**<   Never purge gluons. */ 
    beforeandafter = 1,  /**<   Purge gluons before and after the cascade. */
    everystep = 2        /**<   Purge gluons every step. */
  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  AriadneHandler();

  /**
   * The destructor.
   */
  virtual ~AriadneHandler();
  //@}

public:

  /**
   * The main function to be overwritten by sub-classes. It is called
   * by handle() after storing some information which is then
   * available through simple access functions.
   */
  virtual void cascade();

  /**
   * The AriadneHandler can be used inside the process generation to do
   * so-called CKKW reweighting of the hard sub-process. In this case
   * this function is called after information about the sub-process is
   * made available through the LastXCombInfo base class. Only the
   * function belonging to the primary AriadneHandler for the event to
   * be generated is called. This default implementation of the function
   * simply return one. The current sub-process is mixed together with
   * other processes with a multiplicity of outgoing particles between
   * \a minMult and \a maxMult.
   * @throws CKKWMultiplicityException if the minMult > maxMult or if
   * the number of outgoing particle is outside the specified range.
   */
  virtual double reweightCKKW(int minMult, int maxMult);


  /**
   * Check if the particles in the dipole state passes the cuts from
   * lastXComb.
   */
  bool passCuts(tcDipoleStatePtr state);

public:

  /** @name Functions relating to the running coupling. */
  //@{
  /**
   * Strategy for \f$\alpha_S\f$. 0 means constant a constant
   * \f$\alpha_S\f$ with a value given by alpha0(). 1 means a leading
   * order running \f$\alpha_S\f$ with a \f$\Lambda_{QCD}\f$ given by
   * lambdaQCD(). -1 means take whatever is specified in the current
   * StandardModelBase object.
   */
  inline RunningOption runningCoupling() const {
    return theRunningCoupling;
  }

  /**
   * The constant \f$\alpha_S\f$ to be used if runningCoupling() is 0.
   */
  inline double alpha0() const {
    return theAlpha0;
  }

  /**
   * Return the \f$\alpha_S\f$ to be used for the given \a scale.
   */
  double alphaS(Energy2 scale) const;

  /**
   * The \f$\Lambda_{QCD}\f$ to use in the one loop running
   * \f$\alpha_S\f$ if runningCoupling is 1
   */
  inline Energy lambdaQCD() const {
    return theLambdaQCD;
  }

  /**
   * An internal \f$\alpha_S\f$ object to be used if runningCoupling()
   * is internalRunning.
   */
  inline Ptr<AlphaSBase>::const_pointer internalAlphaS() const {
    return theInternalAlphaS;
  };

  /**
   * Return the flavour thresholds used in calculating \f$\alpha_S\f$.
   */
  inline const set<Energy2> & flavourThresholds() const {
    return theFlavourThresholds;
  }

  /**
   * The constant \f$\alpha_{EM}\f$ to be used. If zero, use whatever
   * is specified in the current StandardModelBase object.
   */
  inline double alphaEM0() const {
    return theAlphaEM0;
  }

  /**
   * The number of different colour indices available to dipoles.
   */
  inline int nCol() const {
    return theNCol;
  }

  /**
   * The number of possible flavours in a \f$g\to q\bar{q}\f$ splitting.
   */
  inline int nFlav() const {
    return theNFlav;
  }
  //@}

  /**
   * Check if photon emission are switched on or off.
   */
  inline bool photonEmissions() const {
    return thePhotonEmissions;
  }

public:

  /**
   * Calculate the starting scale for the given DipoleState.
   */
  Energy startingScale(const DipoleState & state) const;

  /**
   * If DIS-like scattered leptons are found in the given SubProcess,
   * return properly setup corresponding hard remnants. The
   * remnants will belong to the given DipoleState.
   */
  pair<tRemParPtr,tRemParPtr>
  findDISLeptons(SubProcess & sub, DipoleState & state) const;

  /**
   * If findDISLeptons() has found scattered leptons, find also the
   * scattered quarks in the SubProcess and return them as hard
   * remnants. The remnants will belong to the given DipoleState.
   */
  pair<tRemParPtr,tRemParPtr>
  findDISQuarks(pair<tRemParPtr,tRemParPtr> leptons,
		SubProcess & sub, DipoleState & state) const;

  /**
   * Return a list of intermediate s-channel resonances in the given
   * SubProcess. Resonances which have decayed into further resonances
   * will come before their children in the list.
   */
  tPVector resonances(SubProcess & sub) const;

  /**
   * Return the object responsible for radiating from decay products
   * from coloured resonances.
   */
  tColourResonanceModelPtr colourResonanceModel() const {
    return theColourResonanceModel;
  }

  /**
   * Return the object responsible for radiating from remnants.
   */
  const RemnantModel & remnantModel() const {
    return *theRemnantModel;
  }

  /**
   * Find, create and return the QCD dipoles in the given
   * DipoleState.
   */
  vector<tQCDPtr> findQCDDipoles(DipoleState & state) const;

  /**
   * Find, create and return the electro-magnetic dipoles in the given
   * DipoleState.
   */
  vector<tEMDipPtr> findEMDipoles(DipoleState & state) const;

  /**
   * Check the given DipoleState for consistency.
   */
  bool checkState(DipoleState &, tcEmPtr = tEmPtr());

  /**
   * The cutoff in invariant transverse momentum for QCD emissions.
   */
  inline Energy pTCut() const {
    return thePTCut;
  }

  /**
   * The cutoff in invariant transverse momentum for QED emissions.
   */
  inline Energy pTCutEM() const {
    return thePTCutEM;
  }

  /**
   * The inverse extension of a hadron remnant.
   */
  inline Energy softMu() const {
    return theSoftMu;
  }

  /**
   * The dimension assumed for the extension of a hadron remnant.
   */
  inline double softAlpha() const {
    return theSoftAlpha;
  }

  /**
   * The dimension assumed for the extension of a perturbative remnant.
   */
  inline double hardAlpha() const {
    return theHardAlpha >= 0? theHardAlpha: softAlpha();
  }

  /**
   * The power in the suppression of radiation from extended dipoles.
   */
  inline double beta() const {
    return theBeta;
  }

  /**
   * Get the strategy for purging gluons with too small transverse momentum.
   */
  inline PurgeStrategy purgeStrategy() const {
    return thePurgeStrategy;
  }

  /**
   * Get the factor used to determine how far below the cutoff a gluon may be before it is purged.
   */
  inline double purgeFactor() const {
    return thePurgeFactor;
  }

  /**
   * Return the vector of available reweighting objects.
   */
  inline const vector<DipoleRWPtr> & reweighters() const {
    return theReweighters;
  }

  /**
   * The maximum number of allowed emission in the cascade.
   */
  inline int maxEmissions() const {
    return theMaxEmissions;
  }

  /**
   * Return the vector of available emitter objects.
   */
  inline const vector<EmitterPtr> & emitters() const {
    return theEmitters;
  }

  /**
   * Check if the given dipole state corresponds to a valid
   * Born-level state in a CKKW-L merging and return the scale
   * associated with that state. If not return the negative of the
   * scale associated with the state.
   */
  Energy checkBornState(const DipoleState &) const;

  /**
   * Check if the given dipole state corresponds to a valid state from
   * a matrix element generator in a CKKW-L merging.
   */
  bool checkTreeState(const DipoleState &) const;

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

public:

  /**
   * Exception class used if the outgoing multiplicity is out of range
   * or the maximum multiplicity is larger than the minimum in
   * reweightCKKW.
   */
  struct CKKWMultiplicityException: Exception {};

  /**
   * Exception class used if no reasonable state was reconstructed in
   * the CKKW-L algorithm.
   */
  struct CKKWBornException: Exception {};

  /**
   * Exception class used if the dipole state is not self consistant at
   * any stage of the cascade.
   */
  struct IntegretyException: Exception {};

  /**
   * Exception class indicating no model could be found when requested.
   */
  struct MissingModel: public Exception {};

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
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Return true if this object needs to be initialized before all
   * other objects because it needs to extract PDFs from the event
   * file. This version will return true if runningCoupling() returns
   * internalRunning and no internalAlphaS() object has been given.
   */
  virtual bool preInitialize() const;
  //@}

private:

  /**
   * Function used by the interface.
   */
  Energy minPTCut() const;

  /**
   * Function used by the interface.
   */
  Energy maxLambdaQCD() const;

private:

  /**
   * Strategy for \f$\alpha_S\f$. 0 means constant a constant
   * \f$\alpha_S\f$ with a value given by alpha0(). 1 means a leading
   * order running \f$\alpha_S\f$ with a \f$\Lambda_{QCD}\f$ given by
   * lambdaQCD(). -1 means take whatever is specified in the current
   * StandardModelBase object.
   */
  RunningOption theRunningCoupling;

  /**
   * The constant \f$\alpha_S\f$ to be used if runningCoupling() is
   * noRunning.
   */
  double theAlpha0;

  /**
   * The \f$\Lambda_{QCD}\f$ to use in the one loop running
   * \f$\alpha_S\f$ if runningCoupling is simpleRunning.
   */
  Energy theLambdaQCD;

  /**
   * An internal \f$\alpha_S\f$ object to be used if runningCoupling()
   * is internalRunning.
   */
  Ptr<AlphaSBase>::pointer theInternalAlphaS;

  /**
   * Scale factor used to multiply the emission scales in the argument
   * of \f$\alpha_S\f$.
   */
  double scaleFactor;

  /**
   * The flavour thresholds used in calculating \f$\alpha_S\f$.
   */
  set<Energy2> theFlavourThresholds;

  /**
   * The cutoff in invariant transverse momentum for QCD emissions.
   */
  Energy thePTCut;

  /**
   * The constant \f$\alpha_{EM}\f$ to be used. If zero, use whatever
   * is specified in the current StandardModelBase object.
   */
  double theAlphaEM0;

  /**
   * The cutoff in invariant transverse momentum for QED emissions.
   */
  Energy thePTCutEM;

  /**
   * The number of different colour indices available to dipoles.
   */
  int theNCol;

  /**
   * The number of possible flavours in a \f$g\to q\bar{q}\f$ splitting.
   */
  int theNFlav;

  /**
   * Switch to turn on and off photon emission in the cascade.
   */
  bool thePhotonEmissions;

  /**
   * The inverse extension of a hadron remnant.
   */
  Energy theSoftMu;

  /**
   * The dimension assumed for the extension of a hadron remnant.
   */
  double theSoftAlpha;

  /**
   * The dimension assumed for the extension of a hard remnant.
   */
  double theHardAlpha;

  /**
   * The power in the suppression of radiation from extended dipoles.
   */
  double theBeta;

  /**
   * A list of reweighting objects which may be applied to dipole
   * emissions.
   */
  vector<DipoleRWPtr> theReweighters;

  /**
   * The maximum number of emissions allowed in the cascade.
   */
  int theMaxEmissions;

  /**
   * The vector of EmittorBase objects responsible for generating and
   * performing emissions according to the Dipole Cascade Model and
   * its extentions.
   */
  vector<EmitterPtr> theEmitters;

  /**
   * A vector of BornCheckerBase objects which are used in the CKKW-L
   * algorithm.
   */
  vector<BornCheckerBasePtr> theBornCheckers;

  /**
   * The object responsible for choosing the starting scale of the
   * shower.
   */
  ScaleSetterPtr theScaleSetter;

  /**
   * The object responsible for identifying DIS-like event.
   */
  DISFinderPtr theDISFinder;

  /**
   * The object responsible for identifying DIS-like event.
   */
  ResonanceFinderPtr theResonanceFinder;

  /**
   * The Object responsible for identifying QCD Diploles.
   */
  QCDFinderPtr theQCDFinder;

  /**
   * The Object responsible for identifying electro-magnetic Diploles.
   */
  EMFinderPtr theEMFinder;

  /**
   * The object responsible for radiating from decay products from
   * coloured resonances.
   */
  ColourResonanceModelPtr theColourResonanceModel;

  /**
   * The object responsible for radiating from remnants.
   */
  RemnantModelPtr theRemnantModel;

  /**
   * The object responsible for checking the consistency of a DipoleState.
   */
  CCheckPtr consistency;

  /**
   * If set true, the consistencyh checking of emissions in this
   * dipole state is suspended as the initial state was not
   * consistent.
   */
  bool suspendConsistencyChecks;

  /**
   * The strategy for purging gluons with too small transverse momentum.
   */
  PurgeStrategy thePurgeStrategy;

  /**
   * The factor used to determine how far below the cutoff a gluon may be before it is purged.
   */
  double thePurgeFactor;

  /**
   * A map containing the states generated by reweight CKKW. The map
   * assiciates an XComb to a struct containing the dipole state, the
   * the clustering History and a boolean which specifies if the state
   * has the maximum parton multiplicity or not.
   */
  struct CKKWState {
    HistoryPtr history;
    LorentzRotation rotation;
    bool maxMult;
  };
  map<tXCPtr, CKKWState> theCKKWMap;

private:

  /**
   * Histograms for debugging purposes.
   */
  FactoryBase::tH1DPtr histswing, histall, histglue, histqq, histlam;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AriadneHandler & operator=(const AriadneHandler &);

};

}

#endif /* Ariadne5_AriadneHandler_H */
