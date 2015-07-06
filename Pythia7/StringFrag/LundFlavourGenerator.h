// -*- C++ -*-
#ifndef PYTHIA7_LundFlavourGenerator_H
#define PYTHIA7_LundFlavourGenerator_H
// This is the declaration of the LundFlavourGenerator class.

#include "FragConfig.h"
// #include "LundFlavourGenerator.fh"
// #include "LundFlavourGenerator.xh"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Handlers/FlavourGenerator.h"
#include "ThePEG/Repository/EventGenerator.h"

namespace Pythia7 {

/**
 * The LundFlavouGenerator is the class responsible for the flavour
 * generation according to the Lund scheme of fragmentation. It derives
 * from the FlavourGenerator and overrides the two categories of
 * interfaces : the generateHadron() and the getHadron() methods.<br>
 * 
 * The first one is given one incoming flavour, generates a new quark
 * or diquark, and combines it with the existing flavour, then returns
 * the new hadron and the newly created flavour. The second receives
 * two incoming flavours and tries to combine them into a new hadron.
 * 
 * Different models for flavour production can be selected setting the
 * Switch <a
 * href="LundFlavourGeneratorInterfaces.html#BaryonMode"><code>BaryonMode</code></a>.
 *
 * The <i>popcorn</i> model for baryon production needs extra handling
 * than simple flavour production. It is therefore handled by the
 * LundFlavourHandler. Nevertheless the
 * <code>LundFlavourGenerator</code> encapsulates the methods for the
 * popcorn generation. They are denoted as <code>popXXXX()</code>
 * methods.
 *
 * @see \ref LundFlavourGeneratorInterfaces "The interfaces"
 * defined for LundFlavourGenerator.
 * @see FlavourGenerator
 * @see LundFlavourHandler
 *
 */
class LundFlavourGenerator: public ThePEG::FlavourGenerator {

public:

  /** A vector used for SU6 Weights. */
  typedef vector<double> WeightsVec;

  /** An iterator into the SU6 weights vector. */
  typedef WeightsVec::const_iterator WeightsVecPtr;

  /** A matrix used for SU6 Weights. */
  typedef vector< WeightsVec > WeightsTable;

public:


  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  LundFlavourGenerator();

  /**
   * Copy-constructor.
   */
  LundFlavourGenerator(const LundFlavourGenerator &);

  /**
   * Destructor.
   */
  inline virtual ~LundFlavourGenerator();
  //@}

public:

  /** @name Virtual functions required by the FlavourGenerator class. */
  //@{
  /**
   * Generate a hadron from a quark. Given a quark(antiquark, diquark
   * or antidiquark), choose a quark-antiquark (or
   * antidiquark-diquark) pair. Return (first) a hadron formed by the
   * original quark and the antiquark together with (second) the
   * generated quark. Returns null pointers if the generation failed.
   * @param quark a quark, antiquark, diquark or antidiquark.
   * @return a pair of ParticleData pointers. The \a first is the
   * hadron produced and the \a second is the anti-partner of the
   * (anti-)(di-)quark generated to form the hadron.
   */
  virtual tcPDPair generateHadron(tcPDPtr quark) const;

  /**
   * Get hadron from flavours. Return a hadron with the flavour
   * content given by the (anti-)(di-)quarks in the argument. The
   * arguments are given as ParticleData pointers. The default
   * versions will call the getHadron(long, long).
   * @param inPD1 the first flavour.
   * @param inPD2 the second flavour.
   * @return the corresponding hadron type or null if none could be
   * generated.
   */
  virtual tcPDPtr getHadron(tcPDPtr inPD1, tcPDPtr inPD2) const;

  /**
   * Return a baryon with the flavour content given by the
   * (anti)quarks in the argument.  The arguments are given as
   * particle data pointers. The default versions will call
   * getBaryon(tcPDPtr, tcPDPtr, tcPDPtr). If no corresponding hadron was
   * formed it should return the null pointer.
   * @param q1 the PDG code of the first flavour.
   * @param q2 the PDG code of the second flavour.
   * @param q3 the PDG code of the third flavour.
   * @return the corresponding baryon type or null if none could be
   * generated.
   */
  virtual tcPDPtr getBaryon(long q1, long q2, long q3) const;

  /**
   * Generate a random quark flavour.
   */
  virtual long selectQuark() const;

  /**
   * Generate a random (di)quark flavour.
   */
  virtual long selectFlavour() const;
  //@}

  /**
   * Given an input flavour, \a inPD, selects a method (Q2BaryonDQ(),
   * Q2MesonQ(), PopDQ2MesonDQ() or DQ2BaryonQ()) that will generate
   * the new flavour to be combined with \a inPD to produce the new
   * hadron. \a newPD will be set to the newly created flavour and the
   * hadron will be returned.<br> For the popcorn model the Id number
   * of the diquark curtain-quark is given as \a curtainQid. Therefore
   * \a curtainQid = 0 means either that the popcorn model is not
   * switched on or the final baryon of a popcorn generation has to be
   * produced. The selection of curtain quark is however handled by
   * the LundFlavourHandler.
   */
  virtual PDPtr generateHadron(tcPDPtr inPD, cPDPtr& newPD,
			       long curtainQid=0) const;

public:

  /**
   * Initializes the parameters of the generator that depend either on
   * the the interfaced parameters or on the model chosen for the
   * baryon production.  It calls the setMesonFlavourMixingProbs() and
   * initGenPar() methods.
   */
  void initialize();

  /**
   * Returns the a pointer to a quark (or antiquark) of a qqbar pair 
   * created in the string colour field.
   */
  inline virtual PDPtr getRandomFlavour() const;

  
  /** @name Access parameters and swithces. */
  //@{
  /**
   * Returns the baryon mode.
   */
  inline int BaryonMod() const;

  /**
   * Returns the suppression factor for production of a \f$s\bar{s}\f$
   * pair.
   */
  inline double strangeQSup() const;

  /**
   * Returns true if any selected flavour is extra-suppressed. Currently 
   * returns true if an Eta, Eta' meson is extra-suppressed.
   */
  inline bool extraSuppressed() const;

  //@}

  /** @name Functions dealing with popcorn generation. */
  //@{
  /**
   * When a diquark is first produced in the string fragmentation, if
   * (<code>BaryomMode</code>=2) this method is invoked to know the
   * number of popcorn mesons that should be produced in between the
   * baryon anti-baryon pair.
   */
  inline int PopMesonN(tcPDPtr inDQ) const;

  /**
   * When a diquark is first produced in the string fragmentation, if
   * (<code>BaryomMode</code>=2) this method is invoked to know the
   * number of popcorn mesons that should be produced in between the
   * baryon anti-baryon pair.
   */
  int PopMesonN(long inDQ) const;
 
  /**
   * Given the diquark that initializes the popcorn generation,  
   * selects a curtain quark for the popcorn meson production.
   */
  inline long PopSelectCurtainFlavour(tcPDPtr inDQ) const;

  /**
   * Given the diquark that initializes the popcorn generation,  
   * selects a curtain quark for the popcorn meson production.
   */
  long PopSelectCurtainFlavour(long inDQ) const;

  /**
   * In the process of producing a popcorn meson, the PopDQ2MesonD()
   * function may fail to recombine the curtainQ and the new flavour
   * into a new diquark state. The handler of the popcorn generation
   * is informed about this rejection invoking this method that
   * returns true if the generation failed.
   */
  inline bool PopGenRejected() const;
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

  /**
   * Print out state of an object for debugging purposed.
   */
  void DBprint();

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

private:


  /** @name Meson-production functions. */
  //@{
  /**
   * Is the method for the meson production : \a inQ \f$\rightarrow\f$
   * <code>newQ + newMeson</code>. Receives the incoming quark id, \a
   * inQ, generates a new quark flavour and invokes formMeson() to
   * combine them into a newMeson stored in member variables.
   */
  void Q2MesonQ(long inQ) const;

  /**
   * Given the two constituant quark id numbers, \a inQ and \a newQ,
   * forms a new meson and returns its id.
   */
  long formMeson(long inQ, long newQ) const;

  /**
   * When two flavours combine off to form a meson, this function
   * returns the multiplet of the new meson in function of the heavier
   * quark, \a Hfl, of the two.
   */
  int selectMesonMultiplet(long Hfl) const; 

  /**
   * Given the multiplet of the newly created meson, returns its total
   * angular momentum J.
   */
  inline int MesonSpin(int multipletIdx) const;

  /**
   * Selects the physical state of a light flavour-diagonal meson.
   * Given the flavour (\a Qid) of a qqbar pair and the meson
   * multiplet (\a multipletIdx) of the corresponding meson selects at
   * random the meson state (described internally by an integer number
   * that is returned), according to the different probabilities
   * (theMixingProbabilities) that the \f$q\bar{q}\f$ pair produces
   * the possible diagonal-mesons of the different multiplets.
   */
  inline int RandomFlavourMixing(long Qid, int multipletIdx, double rnd) const;
  //@}

  /** @name Diquark production functions. */
  //@{
  /**
   * When a incoming quark is sent to generateHadron() a choice is
   * made to rather produce a quark or a diquark. This method returns
   * true if a diquark is chosen to be produced.
   */
  inline bool DQproduction() const;

  /**
   * When creating a new diquark (in the Q2BaryonDQ)) this method
   * returns true if the diquark is spin-suppressed.
   */
  inline bool DQspinSup(int DQspin) const;
  //@}

  /** @name Baryaon production functions. */
  //@{
  /**
   * This is a baryon production method for the process \a inQ
   * <code>-> newDQ + newBaryon</code><br> Given a incoming quark, \a
   * inQ generates a new diquark such as the q-diquark combination is
   * acceptable according to the corresponding weight for the baryon
   * production. <b>Warning:</b> Contains a special treatment for
   * diquark produced within the pocorm model.
   */
  void Q2BaryonDQ(long inQ) const;

  /**
   * This is the baryon production method for the process : \a inDQ
   * <code>-> newQ + newBaryon</code><br> Given the incoming diquark,
   * \a inDQ, generates a new flavour such as the weight of the baryon
   * corresponding to the (diquark-quark) combination is acceptable.
   */
  void DQ2BaryonQ(long inDQ) const;
 
  /**
   * Given the quark (\a Q) - diquark(\a DQ) flavour combination and
   * the corresponding SU6 Weights, form the new baryon.
   */
  long formBaryon(long Q, long DQ, WeightsVecPtr theWeights) const;

  /**
   * Given the quark (\a Qid) - diquark (\a DQid) combination, returns
   * from the lookup table <code>theSU6WeightsTable</code> a pointer
   * to the table (internally a vector) holding the relative
   * probabilities for the Q, DQ to join into a baryon of the octet or
   * decuplet multiplet.
   */
  WeightsVecPtr getSU6MultWeightVec(long Qid, long DQid) const;

  /**
   * Given the table selected by the getSU6MultWeightVec() returns the
   * relative probability that the quark-diquark combination forms a
   * baryon in the octet.
   */
  inline double OctetWT(WeightsVecPtr wtIt) const;

  /**
   * Given the table selected by the getSU6MultWeightVec() return the
   * relative probability that the quark-diquark combination forms a
   * baryon in the decuplet.
   */
  inline double DecupletWT(WeightsVecPtr wtIt) const;
  //@}

  /** @name Popcorn functions. */
  //@{
  /**
   * Given a diquark (\a inDQ), returns its weight for the popcorn
   * model in function of its flavour content.
   */
  double PopDQweight(long inDQ) const;

  /**
   * Method of the popcorn meson generation : DQ -> popMeson + DQ'
   * <br> Given the incoming diquark (\a inDQ) and the curtain quark
   * (absCurtainQ), combines the <i>free</i> flavour of inDQ with a
   * new created quark to produce the popcorn Meson -> if the
   * generation is not rejected by the PopConsistentDQ() then
   * reconstructs the outgoing diquark DQ', else a new curtain quark
   * has to be re-selected.<br> <b>Warning:</b> the absolute value of
   * curtainQ (\a absCurtainQ) as input here.
   */
  void PopDQ2MesonDQ(long inDQ, long absCurtainQ ) const;

  /**
   * In popcorn meson production, this method returns false if 
   * the selected curtainQ and the new flavour cannot combine into a 
   * new diquark state. 
   */
  bool PopConsistentDQ(long q1, long q2) const;
  //@}

  /** @name Flavour joining functions. */
  //@{
  /**
   * Performs checks on the two incoming flavours to be combined.
   * Returns false if the flavour combination is not consistent.
   */
  bool consistentJoin(long fl1, long fl2) const;
  //@}

  /** @name Inlined helper functions. */
  //@{
  /**
   * Returns true if the incoming flavour id (\a fl) corresponds to a
   * quark number.
   */
  inline bool isQuark(long fl) const;

  /**
   * Returns true if the incoming flavour id (\a fl) corresponds to a
   * diquark number.
   */
  inline bool isDiquark(long fl) const; 


  /**
   * Given a diquark id number (\a inDQ), returns its heavy (\a hq),
   * light(\a lq) quark and its spin (\a s).
   */
  inline void getDQcontent(int inDQ,  int& hq, int& lq, int& s) const;
  //@}



  /** @name Functions for initialization of parameters. */
  //@{
  /**
   * Computes the suppression factor for Spin1 diquark production.  if
   * the popcorn model is selected scales the quark/diquark
   * suppression factors and sets the diquark weights.<br> <b>Warning:
   * </b>To be called after any change of the interfaced parameters
   * involved.
   */
  void initGenPar();

  /**
   * Initializes the default values of the mixing angles,
   * theMixingAngles, for the 6 meson multiplets.
   */
  void setDefaultMixingAnglesVec();

  /**
   * Given the set of Mixing angles, theMixingAngles for the different
   * meson multiplets (the defaults values or the user defined
   * values), initializes the mixing probabilities theMixingProbVec
   * for the production of diagonal-flavour mesons.<br> <b>Warning:
   * </b>To be called after any change of the mixing angles.
   */
  void setMesonFlavourMixingProbs();

  /**
   * Initializes the table of the non-interfaced SU6 Weights :
   * theSU6WeightsTable.
   */
  void setSU6weights();

  /**
   * Initializes the diquark weights for the popcorn model
   * (DQweight).<br> <b>Warning: </b>Called by initGenPar() after any
   * change of the interfaced parameters.
   */
  void setPopDQweights();
  //@}


private:

  /**
   * The Id number of the generated flavour.
   */
  mutable long newFl;

  /**
   * The Id numbers of the created hadron.
   */
  mutable long theHadron;

  /**
   * True if an Eta, Eta' meson constructed by formMeson()
   * method turns out to be extra-suppressed.
   */
  mutable bool extraSup;

  /**
   * True if the PopDQ2MesonDQ method fails to reconstruct the 
   * outgoing diquark. The generation is rejected and a new curtain quark as 
   * to be selected.   
   */
  mutable bool thePopRejection;

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#BaryonMode"><code>BaryonMode</code></a>
   */
  int theBaryonMod;                       //MSTJ(12)

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#DQsup"><code>DQsup</code></a>
   */
  double DQsup;                            //PARJ(1)

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#S1DQsup"><code>S1DQsup</code></a>
   */
  double IS1DQsup;                         //PARJ(4)

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#SBaryonDecupletSup"><code>SBaryonDecupletSup</code></a>
   */
  double BaryonDecupletSup;                //PARJ(18)  

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#extraEtaSup"><code>extraEtaSup</code></a>
   */
  double extraEtaSup;                      //PARJ(25)

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#extraEtapSup"><code>extraEtapSup</code></a>
   */
  double extraEtapSup;                     //PARJ(26)

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#MesonP"><code>MesonP</code></a>
   */
  double MesonP;                           //PARJ(5) 

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#BssBsup"><code>BssBsup</code></a>
   */
  double BssBsup;                          //PARJ(6)

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#BMsBsup"><code>BMsBsup</code></a>
   */
  double BMsBsup;                          //PARJ(7)

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#sQsup"><code>sQsup</code></a>
   */
  double sQsup;                            //PARJ(2)

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#sDQvQsup"><code>sDQvQsup</code></a>
   */
  double sDQvQsup;                         //PARJ(3)

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#S1lightMesonP"><code>S1lightMesonP</code></a>
   */
  double S1lightMesonP;

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#S1sMesonP"><code>S1sMesonP</code></a>
   */
  double S1sMesonP;

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#S1hqMesonP"><code>S1hqMesonP</code></a>
   */
  double S1hqMesonP;

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#P_S0L1J1"><code>P_S0L1J1</code></a>
   */
  double P_S0L1J1;

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#P_S1L1J0"><code>P_S1L1J0</code></a>
   */
  double P_S1L1J0;

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#P_S1L1J1"><code>P_S1L1J1</code></a>
   */
  double P_S1L1J1;

  /**
   * See documentation of interface <a
   * href="LundFlavourGeneratorInterfaces.html#P_S1L1J2"><code>P_S1L1J2</code></a>
   */
  double P_S1L1J2; 

  /**
   *The mixing angles for the production of flavour-diagonal meson.
   */
  vector<double> theMixingAngles;

  /**
   * The lookup table for diquark weights in the popcorn model, 
   * used to compute the baryon weight when combining the DQ with 
   * a new Q in the Q2BaryonDQ() method. Set by initialize().
   */
  vector<double> DQweight;  

  /**
   * Set in initialization to 3*IS1DQsup.
   */
  double S1DQsup;

  /**
   * Popcorn parameters.
   */
  double POP_sDQvQsup;

  /**
   * Popcorn parameters.
   */
  double POP_S1DQsup;

  /**
   * Popcorn parameters.
   */
  double POP_delta;

  /**
   * Maximun weight for diquark production
   */
  double DQmaxWeight;

  /**
   * Lookup Tables for the SU6 weights.
   */
  WeightsTable theSU6WeightsTable;

  /**
   * Lookup Tables for the mixing probabilities .  
   */
  WeightsTable theMixingProbVec;

  /**
   * Standard Interface
   */
  static ClassDescription<LundFlavourGenerator> initLundFlavourGenerator;

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * Pythia7::LundFlavourGenerator.
 */
template <>
struct BaseClassTrait<Pythia7::LundFlavourGenerator,1>: public ClassTraitsType {
  /** Typedef of the base class of Pythia7::LundFlavourGenerator. */
  typedef FlavourGenerator NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Pythia7::LundFlavourGenerator class and the shared object where it
 * is defined.
 */
template <>
struct ClassTraits<Pythia7::LundFlavourGenerator> 
  : public ClassTraitsBase<Pythia7::LundFlavourGenerator> {
  /** Return the class name.  */
  static string className() { return "Pythia7::LundFlavourGenerator"; }
  /**
   * Return the name of the shared library to be loaded to get access
   * to the Pythia7::LundFlavourGenerator class and every other class
   * it uses (except the base class).
   */
  static string library() { return "libP7String.so"; }
};

/** @endcond */

}

#include "LundFlavourGenerator.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "LundFlavourGenerator.tcc"
#endif

#endif /* PYTHIA7_LundFlavourGenerator_H */




