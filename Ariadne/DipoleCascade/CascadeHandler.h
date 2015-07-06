// -*- C++ -*-
#ifndef ARIADNE_CascadeHandler_H
#define ARIADNE_CascadeHandler_H
//
// This is the declaration of the CascadeHandler class.
//

#include "Ariadne/Config/Ariadne.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/Utilities/Triplet.h"
#include "ThePEG/StandardModel/AlphaSBase.h"
#include "CascadeHandler.fh"
#include "MECorrBase.fh"
#include "ReweightBase.fh"
#include "DipoleState.fh"

namespace Ariadne {

using namespace ThePEG;

/**
 * The CascadeHandler class inherits form the ThePEG::CascadeHandler
 * base class and implements the Dipole Cascade Model for partonic
 * cascades.
 *
 * @see \ref CascadeHandlerInterfaces "The interfaces"
 * defined for CascadeHandler.
 */
class CascadeHandler: public ThePEG::CascadeHandler {

public:

  /**
   * A vector of matrix element correction objects.
   */
  typedef vector<MECPtr> MECVector;

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

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  CascadeHandler();

  /**
   * The copy constructor.
   */
  CascadeHandler(const CascadeHandler &);

  /**
   * The destructor.
   */
  virtual ~CascadeHandler();
  //@}

public:

  /**
   * The main function to be overwritten by sub-classes. It is called
   * by handle() after storing some information which is then
   * available through simple access functions.
   */
  virtual void cascade();

  /**
   * The CascadeHandler can be used inside the process generation to do
   * so-called CKKW reweighting of the hard sub-process. In this case
   * this function is called after information about the sub-process is
   * made available through the LastXCombInfo base class. Only the
   * function belonging to the primary CascadeHandler for the event to
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
  inline RunningOption runningCoupling() const;

  /**
   * The constant \f$\alpha_S\f$ to be used if runningCoupling() is 0.
   */
  inline double alpha0() const;

  /**
   * The \f$\Lambda_{QCD}\f$ to use in the one loop running
   * \f$\alpha_S\f$ if runningCoupling is 1
   */
  inline Energy lambdaQCD() const;

  /**
   * An internal \f$\alpha_S\f$ object to be used if runningCoupling()
   * is internalRunning.
   */
  inline Ptr<AlphaSBase>::const_pointer internalAlphaS() const;

  /**
   * The constant \f$\alpha_{EM}\f$ to be used. If zero, use whatever
   * is specified in the current StandardModelBase object.
   */
  inline double alphaEM0() const;

  /**
   * The number of different colour indices available to dipoles.
   */
  inline int nCol() const;

  /**
   * The number of possible flavours in a \f$g\to q\bar{q}\f$ splitting.
   */
  inline int nFlav() const;
  //@}

  /**
   * Check if photon emission are switched on or off.
   */
  inline bool photonEmissions() const;

public:

  /**
   * The cutoff in invariant transverse momentum for QCD emissions.
   */
  inline Energy pTCut() const;

  /**
   * The cutoff in invariant transverse momentum for QED emissions.
   */
  inline Energy pTCutEM() const;

  /**
   * The inverse extension of a hadron remnant.
   */
  inline Energy softMu() const;

  /**
   * The dimension assumed for the extension of a hadron remnant.
   */
  inline double softAlpha() const;

  /**
   * The dimension assumed for the extension of a perturbative remnant.
   */
  inline double hardAlpha() const;

  /**
   * The power in the suppression of radiation from extended dipoles.
   */
  inline double beta() const;

  /**
   * Return the vector of available reweighting objects.
   */
  inline const vector<DipoleRWPtr> & reweighters() const {
    return theReweighters;
  }

  /**
   * Return the vector of available matrix element reweighting objects.
   */
  inline const MECVector & MECorrectors() const;

  /**
   * The maximum number of allowed emission in the cascade.
   */
  inline int maxEmissions() const;

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

  /**
   * Exception class used if the outgoing multiplicity is out of range
   * or the maximum multiplicity is larger than the minimum in
   * reweightCKKW.
   */
  struct CKKWMultiplicityException: Exception {};

  /**
   * Exception class used if the dipole state is not self consistant at
   * any stage of the cascade.
   */
  struct IntegretyException: Exception {};

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

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
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointers to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();

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
   * A list of leading-order matrix-element correction objects which
   * may be applied to dipole emissions.
   */
  MECVector theMECorrectors;

  /**
   * The maximum number of emissions allowed in the cascade.
   */
  int theMaxEmissions;

  /**
   * A map containing the states generated by reweight CKKW. The map
   * assiciates an XComb to a struct containing the dipole state, the
   * rotation, the last pt2 and a boolean which specifies if the state
   * has the maximum parton multiplicity or not.
   */
  struct CKKWState {
    DipoleStatePtr state;
    LorentzRotation rotation;
    Energy2 pt2;
    bool maxMult;
  };
  map<tXCPtr, CKKWState> theCKKWMap;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<CascadeHandler> initCascadeHandler;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CascadeHandler & operator=(const CascadeHandler &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of CascadeHandler. */
template <>
struct BaseClassTrait<Ariadne::CascadeHandler,1> {
  /** Typedef of the first base class of CascadeHandler. */
  typedef CascadeHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the CascadeHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<Ariadne::CascadeHandler>
  : public ClassTraitsBase<Ariadne::CascadeHandler> {
  /** Return a platform-independent class name */
  static string className() { return "Ariadne::CascadeHandler"; }
  /** Return the name of the shared library be loaded to get
   *  access to the CascadeHandler class and every other class it uses
   *  (except the base class). */
  static string library() { return "libArCascade.so"; }
};

/** @endcond */

}

#include "CascadeHandler.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CascadeHandler.tcc"
#endif

#endif /* ARIADNE_CascadeHandler_H */
