// -*- C++ -*-
#ifndef PYTHIA7_LundFragHandler_H
#define PYTHIA7_LundFragHandler_H
// This is the declaration of the LundFragHandler class.

#include "FragConfig.h"
#include "ThePEG/Handlers/HadronizationHandler.h"
#include "LundFragHandler.xh"
#include "LundFlavourGenerator.h"
#include "Oriented.h"
#include "OrientedIndex.h"
#include "ThePEG/Handlers/PtGenerator.h"
#include "ThePEG/Handlers/ZGenerator.h"
#include "ThePEG/Handlers/FlavourGenerator.h"
#include "ThePEG/Handlers/ClusterCollapser.h"
#include "EndPoint.h"
#include "Hadron.h"


namespace Pythia7 {

/**
 * The LundFragHandler is the main class of the string fragmentation
 * modules. It is responsible for handling the hadronization phase
 * according to the Lund fragmentation scheme. It inherits from
 * HadronizationHandler the methods to communicate with the
 * EventHandler and the EventGenerator.  It also derives from Oriented
 * to fix the orientation taken for each step in the fragmentation
 * procedure.
 *
 * Its tasks is : to look through the <code>tagged</code> particles
 * present in the current Step, to extract the Strings to be
 * hadronized into a list of Particles. At the end of the procedure,
 * it asks the EventHandler for a copy of the current Step
 * to add the newly created hadrons.
 *  
 * The LundFragHandler is responsible for creating the current String
 * object to fragment. The fragmentation algorithm is implemented as a
 * step-by-step updating procedure of EndPoint.  The LundFragHandler
 * contains three <code>EndPoint</code>s: two to describe the last
 * (right, left) end-point left over in the previous steps, and the
 * current end-point created during the current step.
 *
 * To achieve its job the LundFragHandler makes use of the
 * LundPtGenerator, LundZGenerator and LundFlavourGenerator classes.
 *
 * @see \ref LundFragHandlerInterfaces "The interfaces"
 * defined for LundFragHandler. 
 * @see HadronizationHandler
 * @see Oriented
 * @see String
 * @see EndPoint
 */
class LundFragHandler:public HadronizationHandler, public Oriented {

public:

  /** Shorthand alias. */
  typedef OrientedIndex OIndex;

  /** A pair of energy fractions. */
  typedef pair<double, double> Xhat;
  /** A vector of energy fraction pairs. */
  typedef vector<Xhat> XhatVector;
  /** A list of Hadron object. */
  typedef list<Hadron> Buffer;
  /** Iterator into a list of Hadron object. */
  typedef Buffer::iterator BufferIt;


public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  LundFragHandler();

  /**
   * Copy-constructor.
   */
  LundFragHandler(const LundFragHandler &);

  /**
   * Destructor.
   */
  ~LundFragHandler();
  //@}

  /** @name Virtual functions to be implemented by concrete sub-classes. */
  //@{
  /**
    * The main function called by the EventHandler class to
    * perform a step.  Extracts the color singlet strings in the
    * tagged particles vector, and send them to the Hadronize method
    * to be fragmented.  When completed get a copy of the current Step
    * from the EventHandler and insert the newly created
    * particles to form a new step in the EventRecord.
    * @param eh the EventHandler in charge of the Event generation.
    * @param tagged if not empty these are the only particles which should
    * be considered by the StepHandler.
    * @param hint a Hint object with possible information from previously
    * performed steps.
    * @throws Veto if the StepHandler requires the current step to be
    * discarded.
    * @throws Stop if the generation of the current Event should be stopped
    * after this call.
    * @throws Exception if something goes wrong.
    */
  virtual void handle(EventHandler & eh, const tPVector & tagged,
		      const Hint & hint);
  //@}

  /**
   * The major method in the string fragmentation administration. 
   * Given a vector of Particles that form the string, it returns the list 
   * of newly produced particles. Its task is divided in three main parts : <br>
   * 1) Set-up the string calling the initHadronization method<br>
   * 2) call getHadron() method to produce a new hadron while the 
   *    remaining energy of the string is above a minimum threshold, 
   *    otherwise call the finalTwoHadron() method the produce 
   *    the  last 2 hadron in the procedure.
   */
  ParticleList Hadronize(const tcPVector& );


  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();


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

public: 

  /** @name Access helper generator objects. */
  //@{
  /**
   * Set pointer to LundPtGenerator.
   */
  inline void PtGen(PtGeneratorPtr);

  /**
   * Get pointer to LundPtGenerator.
   */
  inline PtGeneratorPtr PtGen() const;

  /**
   * Set pointer to LundFlavourGenerator.
   */
  inline void FlavourGen(FlavourGeneratorPtr);

  /**
   * Get pointer to LundFlavourGenerator.
   */
  inline FlavourGeneratorPtr FlavourGen() const;

  /**
   * Set pointer to LundZGenerator.
   */
  inline void ZGen(ZGeneratorPtr);

  /**
   * Get pointer to LundZGenerator.
   */
  inline ZGeneratorPtr ZGen() const;
  //@}

  /** @name Access interfaced parameters. */
  //@{
  /**
   * Get maximum number of attempts.
   */
  inline long maxLoop() const;

  /**
   * Set maximum number of attempts.
   */
  inline void maxLoop(long n);

  /**
   * See documentation of interface <a
   * href="LundFragHandlerInterfaces.html#Wmin0"><code>Wmin0</code></a>
   */
  inline Energy Wmin0() const;

  /**
   * See documentation of interface <a
   * href="LundFragHandlerInterfaces.html#k"><code>k</code></a>
   */
  inline double k() const;

  /**
   * See documentation of interface <a
   * href="LundFragHandlerInterfaces.html#delta"><code>delta</code></a>
   */
  inline double Delta() const;

  /**
   * See documentation of interface <a
   * href="LundFragHandlerInterfaces.html#m_0"><code>m_0</code></a>
   */
  inline Energy m0() const;

  /**
   * See documentation of interface <a
   * href="LundFragHandlerInterfaces.html#d_0"><code>d_0</code></a>
   */
  inline InvEnergy4 d0() const;

  /**
   * Returns the suppression factor for production of a \f$s\bar{s}\f$
   * pair (from the flavour generator).
   */
  inline double SqRatio() const;
  //@}

  
  /** @name Print functions for debugging. */
  //@{
  /**
   * Print end points.
   */
  void showEP() const;

  /**
   * Print Xhat values for a string region.
   */
  void echoXhat(cStringRegionPtr );
  //@}

protected:


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
  inline virtual void doinit() throw(InitException);

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  //@}

protected: 

  /** @name Main internal functions for the hadronization. */
  //@{
  /**
   * Set the LundFraghandler members to their defaults value
   */
  void resetHandler();

  /**
   * Find the string type (open, closed) of the incoming vector of
   * particles
   */
  void initHadronization(const tcPVector& );

  /**
   * Copy partons which should be fragmented and recombine those which
   * are too close together.
   */
  PVector copyRecombine(tcPPtr first, const tcPVector & inPVec);

  /**
   * Given the vector of particles, find the correct initial particle
   * to send to the String constructor to create the open String and
   * initialize the first <code>EndPoint</code>s
   */
  void initOpenString(const tcPVector& );  

  /**
   * Given the vector of particles, find the correct initial particle
   * to send to the String constructor to create the closed String and
   * initialize the first <code>EndPoint</code>s
   */
  void initClosedString(const tcPVector& );  


  /**
   * Perform a step in the fragmentation iterative process, resulting 
   * by the creation of a new hadron. This is the method that fixes 
   * the step's side, that will be used by all oriented classes.
   */
  void getHadron();

  /**
   * Updates the CurrentString variables after a step has been preformed 
   * (i.e. a new break-point created and a new hadron produced) before  
   * starting a new step.  
   */
  void loopBack();

  /**
   * Produce the last two hadron in the string fragmentation scheme.
   * Invoked when the String remaining energy is below Wmin. At that 
   * point the iterative procedure is stopped andd the last two hadrons
   * are produced 
   */
  void finalTwoHadrons();

  /**
   * Main method responsible for finding the string region where a 
   * solution for the new generated breakpoint can be found.
   */
  void Stepping();

  /**
   * Invoked by Stepping() to solve the (Gamma, m2) equation system in
   * the current StringRegion reached by the Stepping procedure.
   */
  void solveGammaM2System();

  /**
   * Handle the steps along the StringRegion map.
   */
  void Step();
  /**
   * Handle the steps along the StringRegion map.
   */
  void stepDown();

  /**
   * Calculate the momentum vector for two string regions.
   */
  LorentzMomentum p0(cStringRegionPtr first, cStringRegionPtr last);
  

  /**
   * Used by <code>initClosedString</code> to break a closed string
   * and find the rightmost, leftmost end-points, for the resulting 
   * open-like string. 
   */
  cPPtr selectBreakup(const tcPVector&);

  /**
   * Used by <code>initClosedString</code> to break a closed string
   * and find the rightmost, leftmost end-points, for the resulting 
   * open-like string. 
   */
  void pickFirstEPts();

  //@}

  /** @name Functions for creating the final two hadrons. */
  //@{
  /**
   * Used by the finalTwoHadrons() to find the final region when the 
   * two last End-point are not in the same string-region
   */
  void setupCommonFinalRegion();

  /**
   * Used by finalTwoHadrons() to compute the kinematics of the two 
   * final hadrons 
   */
  void solveKinematics(); 

  /**
   * Return the final string-region where the second end-point is produced.
   */
  inline cStringRegionPtr finalSR() const;

  /**
   * Return the first end-point. 
   */
  inline  const EndPoint& firstEP() const;

  /**
   * Return the second end-point. 
   */
  inline  EndPoint& secondEP();

  /**
   * Get Ref to the second Hadron.
   */
  inline  Hadron& secondH();
  //@}



  /** @name Helper and access functions. */
  //@{
  /**
   * Given the masses of the last two end-point and the current one ,
   * compute the minimum energy above which a \f$q\bar{q}\f$ pair can
   * be produced.
   */
  Energy2 Wmin2() const;

  /**
   * Return true if the string remaining energy is large enough 
   * to produce a new \f$q\bar{q}\f$ pair. 
   */
  inline bool enoughE() const; 

  /**
   * Return true if the current string is a simple \f$q\bar{q}\f$
   * string.
   */
  inline bool AqqbarSystem() const;

  /**
   * Used by the Step() method. Return true if a step taken on the 
   * string region map ends-up in an inconsistent string-region, 
   * that will lead to uncrossed the left-right sequences.
   */
  inline bool inconsistentBreakupRegions() const;

  /**
   * Return the number of string-regions of the currentString.
   */
  inline int  nSR() const;

  /**
   * EndPoints accessor.
   */
  inline const EndPoint& lastEP() const;

  /**
   * EndPoints accessor. 
   */
  inline const EndPoint& lastOppEP() const;

  /**
   * EndPoints accessor. 
   */
  inline EndPoint& getLastEP();

  /**
   * EndPoints accessor. 
   */
  inline EndPoint& getLastOppEP();

  /**
   * Return the current string
   */
  inline cStringRegionPtr CurrentSR() const;

  /**
   * Access to the forward momentum fraction available in the
   * current StringRegion.
   */
  inline double CurrentXremf() const;

  /**
   * Access to the backward momentum fraction available in the current
   * StringRegion.
   */
  inline double CurrentXremb() const;

  /**
   * Get ref to the forward Xhat fractions given the index (int) of
   * the string-region axis. (cf. StrinRegion).
   */
  inline double& Xhatfwd(int );

  /**
   * Get ref to the backward Xhat fractions given the index (int) of
   * the string-region axis. (cf. StrinRegion).
   */
  inline double& Xhatbwd(int );

  /**
   * Set the forward and backward Xhat fractions 
   * in the current StrinRegion.
   */
  inline void setXhat(double, double);

  /**
   * Return true when a solution has been found for the position 
   * of the new breakup.
   */
  inline bool aSolution() const;
  //@}

  /** @name Functions for generating and storing hadrons. */
  //@{
  /**
   * Store a hadron.
   */
  void store(Hadron& , int Dir=Oriented::Dir());

  /**
   * Clear hadron storage.
   */
  inline void clearBuffer();

  /**
   * Create a list of stored hadrons.
   */
  ParticleList createParticleList();

  /**
   * Calls LundFlavourGenerator::generateHadron(tcPDPtr, cPDPtr&, long).
   */
  virtual tcPDPtr generateHadron(tcPDPtr inPDPtr, cPDPtr& newPDPtr,
				 long curtainQid=0);

  /**
   * Calls LundFlavourGenerator::getHadron(tcPDPtr, tcPDPtr).
   */
  virtual tcPDPtr getHadron(tcPDPtr inPD1, tcPDPtr inPD2);
  //@}


private : 

  /**
   * The LundPtGenerator.
   */
  PtGeneratorPtr thePtGen;

  /**
   * The LundFlavourGenerator.
   */
  FlavourGeneratorPtr theFlGen;  

  /**
   * The LundZGenerator.
   */
  ZGeneratorPtr theZGen;

  /**
   * The object used to avoid too small strings in the hadronization.
   */
  ClusterCollapserPtr theCollapser;

  /**
   * See documentation of interface <a
   * href="LundFragHandlerInterfaces.html#Wmin0"><code>Wmin0</code></a>
   */
  Energy pWmin0;          // PARJ(33) - PARJ(34) 

  /**
   * See documentation of interface <a
   * href="LundFragHandlerInterfaces.html#k"><code>k</code></a>
   */
  double pK;              // PARJ(36)

  /**
   * See documentation of interface <a
   * href="LundFragHandlerInterfaces.html#delta"><code>delta</code></a>
   */
  double pDelta;          // PARJ(37)

  /**
   * See documentation of interface <a
   * href="LundFragHandlerInterfaces.html#m_0"><code>m_0</code></a>
   */
  Energy pM0;             // PARJ(32)


  /**
   * The effective cut-off in squared mass, below which partons may be
   * recombined to simplify (machine precision limited) kinematics of
   * string fragmentation. (Default chosen to be of the order of a
   * light quark mass, or half a typical light meson mass.)
   */
  Energy2 m2min;         // PARU(12)

  /**
   * The effective cut-off in squared mass, below which partons may be
   * recombined to simplify (machine precision limited) kinematics of
   * string fragmentation. (Default chosen to be of the order of a
   * light quark mass, or half a typical light meson mass.) m2mini is
   * copied from m2min for each string to be hadronized and will be
   * increased if the hadronization fails too often.
   */
  Energy2 m2mini;

  /**
   * The effective angular cut-off in radians for recombination of
   * partons, used in conjunction with m2min.
   */
  double angmin;         // PARU(13)

  /**
   * The effective angular cut-off in radians for recombination of
   * partons, used in conjunction with m2min. angmini is
   * copied from angmin for each string to be hadronized and will be
   * increased if the hadronization fails too often.
   */
  double angmini;

  /**
   * Used to parametrize the probability for reverse rapidity ordering 
   * of the final two hadrons.
   */
  InvEnergy4 pd0;         // PAR(38)

  /**
   * Defines the maximun mumber of allowed attemps to fragment the CurrentString
   */
  long MaxLoop;


  /**
   * The current string to fragment
   */
  StringPtr theCurrentString;

  /**
   * The right end-point.
   */
  EndPoint  theRightEP;

  /**
   * The left end-point. 
   */
  EndPoint  theLeftEP;

  /**
   * The current end-points 
   */
  EndPoint  CurrentEP;

  /**
   * The final string-region in the final 2 hadrons procedure. 
   */
  cStringRegionPtr thefinalSR;         

  /**
   * Newly created ParticleData type.
   */
  cPDPtr newCreatedPD;        

  /**
   * Newly created Hadrons.
   */
  Hadron newHadron;

  /**
   * Newly created Hadrons.
   */
  Hadron lastHadron;


  /**
   * Vector of forward Xhat coordinates
   */
  XhatVector XhatFwdVector;

  /**
   * Vector of backward Xhat coordinates
   */
  XhatVector XhatBwdVector;

  /**
   * Used to keep track of total momentum available in breakups.
   */
  LorentzMomentum Pzero;

  /**
   * True as long as there is enough energy left in the string.
   */
  bool Estatus;

  /**
   * True as long as there is a solution to the gamma equation.
   */
  bool GammaM2Solution;

  /**
   * The number of attempts so far to fragment the current string.
   */
  long ntry;

  /**
   * Closed string first breakup 
   */
  cPPtr breakup;

  /**
   * Hadron storage.
   */
  Buffer theBuffer;

  /**
   * Hadron storage iterator.
   */
  BufferIt currentBufferIt; 

  /**
   * The rotation to boost back ftom the current string rest frame.
   */
  LorentzRotation cmr;

  /**
   * Interface description
   */
  static ClassDescription<LundFragHandler> initLundFragHandler;

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * Pythia7::LundFragHandler.
 */
template <>
struct BaseClassTrait<Pythia7::LundFragHandler,1>: public ClassTraitsType {
  /** Typedef of the base class of Pythia7::LundFragHandler. */
  typedef HadronizationHandler NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Pythia7::LundFragHandler class and the shared object where it
 * is defined.
 */
template <>
struct ClassTraits<Pythia7::LundFragHandler> 
  : public ClassTraitsBase<Pythia7::LundFragHandler> {
  /** Return the class name.  */
  static string className() { return "Pythia7::LundFragHandler"; }
  /**
   * Return the name of the shared library to be loaded to get access
   * to the Pythia7::LundFragHandler class and every other class
   * it uses (except the base class).
   */
  static string library() { return "libP7String.so"; }
};

/** @endcond */

}

#include "LundFragHandler.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "LundFragHandler.tcc"
#endif

#endif /* PYTHIA7_LundFragHandler_H */





