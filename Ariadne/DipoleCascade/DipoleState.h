// -*- C++ -*-
#ifndef ARIADNE_DipoleState_H
#define ARIADNE_DipoleState_H
//
// This is the declaration of the DipoleState class.
//

#include "CascadeBase.h"
#include "DipoleState.fh"
#include "Emitter.h"
#include "SoftRemnant.h"
#include "String.h"
#include "HardSubSys.h"
#include "ThePEG/Utilities/ObjectIndexer.h"

namespace Ariadne {

/**
 * DipoleState contains all Ariadne::Parton, Ariadne::Dipole and
 * Ariadne::String objects needed to describe a partonic state within
 * the Dipole Cascade Model.
 */
class DipoleState: public CascadeBase {

public:

  /**
   * A set of CascadeBase pointers.
   */
  typedef list<CascadeBasePtr> BaseSet;

  /**
   * A set of Emitter pointers.
   */
  typedef list<tEmiPtr> EmitterSet;

  /**
   * A set of String pointers.
   */
  typedef list<tStrPtr> StringSet;

  /**
   * A pair of partons.
   */
  typedef pair<tSoftRemPtr,tSoftRemPtr> tRemPair;

  /**
   * A Selector for DipoleStates
   */
  typedef Selector<DipoleStatePtr> DipoleStateSelector;

  /**
   * Enumerate the different types of sub-processes.
   */
  enum ProcessType {
    unknown,      /**< Used if not cascading a primary sub-process. */
    annihilation, /**< Lepton annihilation. */
    leptonlepton, /**< Virtual photon-photon scattering. */
    leptonhadron, /**< Virtula photon-hadron scattering. */
    hadronlepton, /**< Hadron-virtual photon scattering. */
    hadronhadron, /**< Virtual photon-photon scattering. */
    DIPSYlephad,  /**< Virtula photon-hadron scattering. */
    DIPSYhadlep,  /**< Hadron-virtual photon scattering. */
    DIPSYhadhad   /**< Hadron-hadron collision. */
  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline DipoleState(tHandlerPtr hdl = tHandlerPtr());

  /**
   * The copy constructor.
   */
  inline DipoleState(const DipoleState &);

  /**
   * The destructor.
   */
  virtual ~DipoleState();
  //@}

public:

  /** @name The main functions controlling the cascade. */
  //@{
  /**
   * Initialize this state. Using information about the hard
   * sub-process, initialize the DipoleState, preparing for the
   * cascade.
   *
   * @param sub the SubProcess
   *
   * @return the Lorentz rotation used to boost to the hadronic
   * center-of-mass system.
   */
  LorentzRotation init(tSubProPtr sub);

  /**
   * Init simple final-state cascade. Using the outgoing partons in \a
   * out to initialize a final-state cascade.
   *
   * @return false if no cascade was possible.
   */
  bool init(const tPVector & out);

  /**
   * Select an Emitter. Go through the list of Emitter objects and
   * tell them to calculate the invariant transverse momentum of the
   * next emission and select the Emitter with largest transverse
   * momentum. Subsequent calls to selected() will return the selected
   * Emitter.
   *
   * @param pt2min the minimum allowed transverse momentum.
   *
   * @param pt2max the maximum allowed transverse momentum.
   *
   * @return the generated invariant squared transverse momentum. If
   * ptmin is returned, no emission was generated and subsequent calls
   * to selected() will return null.
   */
  Energy2 select(Energy2 pt2min, Energy2 pt2max);

  /**
   * Perform an emission.
   *
   * @return false if no emission was possible.
   */
  bool perform();

  /**
   * Output the result of a simple final-state cascade. The new
   * partons are created, boosted and inserted in the Step provided.
   *
   * @param inital the initial list of partons before the cascade.
   *
   * @param rot the LorentzRotation by which the final partons will be
   * bosted.
   *
   * @param step the step in which to insert the produced partons.
   */
  void fill(const tPVector & initial, const LorentzRotation & rot, tStepPtr step);

  /**
   * Output the result of a cascade starting from a full sub
   * process. The resulting partons are created, boosted and inserted
   * in the Step provided.
   *
   * @param sub the initial SubProcess object.
   *
   * @param rot the LorentzRotation by which the final partons will be
   * bosted.
   *
   * @param step the step in which to insert the produced partons.
   */
  void fill(tSubProPtr sub, const LorentzRotation & rot, tStepPtr step);
  //@}

  /** @name Simple access functions. */
  //@{
  /**
   * The Emitter selected to perform the next emission.
   */
  inline tEmiPtr selected() const;  

  /**
   * Set the Emitter selected to perform the next emission.
   */
  inline void selected(const tEmiPtr &);  

  /**
   * Get the pair of remnant partons, if any.
   */
  inline const tRemPair & remnants() const;

  /**
   * Set the pair of remnant partons, if any.
   */
  inline void remnants(const tRemPair &);

  /**
   * Replace an \a old remnant parton with a new \a p.
   * @return false if \a old was not an remnant parton.
   */
  bool replaceRemnant(tSoftRemPtr old, tSoftRemPtr p);

  /**
   * Access information abot the current hard sub-system.
   */
  inline const HardSubSys & hardSubSys() const;

  /**
   * Access information abot the current hard sub-system.
   */
  inline HardSubSys & hardSubSys();

  /**
   * Return the squared invariant mass of the hadronic center-of-mass
   * system.
   */
  Energy2 sTot() const;

  /**
   * The Process type.
   */
  inline ProcessType pType() const;

  /**
   * The number of emissions made in this state.
   */
  inline int ne() const;

  /**
   * Return the inverse extension \f$\mu\f$ of the remnant from
   * particle \a p.
   */
  Energy getSoftMu(tcPPtr p) const;

  /**
   * Return the inverse extension \f$\mu\f$ of a hard probe \a p.
   */
  Energy getHardMu(tcPPtr p) const;

  /**
   * Return the dimensionality \f$\alpha\f$ of the extension of
   * remnant from particle \a p.
   */
  double getSoftAlpha(tcPPtr p) const;

  /**
   * Return the dimensionality of the extension of a hard probe \a p.
   */
  double getHardAlpha(tcPPtr p) const;

  /**
   * Return the extendedness of the remnant particles.
   */
  inline const pair<Energy,Energy> & mus() const;

  /**
   * The incoming particles giving the hadronic system.
   */
  inline tPPair particles() const;

  /**
   * The scattered lepton(s) in case this is a DIS event.
   */
  inline tPPair scatteredLeptons() const;

  /**
   * The scattered quark(s) in case this is a DIS event.
   */
  inline tPPair scatteredQuarks() const;
  //@}

  /** @name Functions relating to the book keeping of included objects. */
  //@{
  /**
   * Create a new object to be included in this state.
   */
  template <typename Class>
  inline typename Ptr<Class>::pointer create();

  /**
   * Create a new Parton of the given type to be included in this state.
   */
  ParPtr create(PDPtr type);

  /**
   * Create a new Parton from the given Particle to be included in
   * this state.
   */
  ParPtr create(tcPPtr p);

  /**
   * Remove the given parton. Note that it is only removed from the
   * list of partons. It is still remembered as an object belonging to
   * this state.
   */
  void remove(tParPtr p);

  /**
   * Remove the given emitter. Note that it is only removed from the
   * list of emitters. It is still remembered as an object belonging to
   * this state.
   */
  void removeEmitter(tEmiPtr p);

  /**
   * Remove the given string. Note that it is only removed from the
   * list of strings. It is still remembered as an object belonging to
   * this state.
   */
  void removeString(tStrPtr p);

  /**
   * Clone this state. Also cloning all included objects and fixing up
   * their inter-dependence.
   */
  DipoleStatePtr fullclone() const;

  /**
   * Clone this state. Also cloning all included objects and fixing up
   * their inter-dependence. After the call, \a trans will contain the
   * translation map between the objects in the old and new
   * DipoleState.
   */
  DipoleStatePtr fullclone(TranslationMap & trans) const;

  /**
   * Clone this state. Also cloning all included objects but do not
   * fix up their inter-dependence, ie. the cloned objects in the new
   * DipoleState may still point to objects in the old state. Before
   * the new state can be used, its postclone() function must be
   * called with the \a trans object as argument. The \a trans object
   * will contain the translation map between the objects in the old
   * and new DipoleState.
   */
  DipoleStatePtr preclone(TranslationMap & trans) const;

  /**
   * Fix up the inter-dependence among the included objects in this
   * DipoleState which was produced by the preclone() function. The
   * same \a trans object which was used in the preclone() call must
   * be given here as argument.
   */
  void postclone(const TranslationMap & trans) const;
  //@}

public:

  /**
   * Construct all possible histrories to this dipole state
   * and select one radomly according to the product of the
   * splitting functions. A PT ordered history is selected if
   * possible The selected state is saved in theSlectedHistory. 
   * theSelected will be set to the modified dipole. Returns false 
   * if no history was found.
   */
  bool constructHistory(int steps);

  /**
   * Returns the product of coupling constants for at the constructed scales.
   */
  double couplingProduct(int steps);

  /**
   * Returns the product of the pdf ratios that would have been used in
   * the undone emissions.
   */
  double PDFRatioProduct(int steps);

  /** 
   * Calculates whether or not an event should be vetoed according to
   * the Sudakov Veto algorithm. A return value of true means veto the
   * event. This function does not take in to account the veto that may
   * occur when emissions from the state given by the maxtrix element is
   * above the matrix element cutoff.  This function destroys the
   * information in theSelectedHistory. Any other reweighing should be
   * done before calling this function.
   */
  bool sudakovVeto(int steps);

  /**
   * Find a matrix-element correction object for the given \a dipole;
   */
  tcMECPtr findMECorr(tcEmiPtr dipole) const;

  /**
   * Calculate the total momentum of the dipole state.
   */
  LorentzMomentum totalMomentum() const;

  /**
   * Return the constructed pt2 of the last emission.
   */
  inline Energy2 constructedPT2() const;

protected:

  /**
   * A function that constructs all possible histories.
   */
  void allHistories(int steps);

  /**
   * Returns true if the exists at least one ordered history.
   */
  bool hasOrderedHistory(int steps);

  /**
   * Removes all unordered histories. Returns true if there are
   * any states left.
   */
  bool removeUnorderedHistories(int steps);

  /**
   * Select one history among the ones constructed. Removes all of 
   * the states in theHistory and saves the selected one in 
   * theSelectedHistory. Returns false if no possible history 
   * was found. Modifies the pT2 values to make the cascade history
   * ordered if it is unordered.
   */
  bool selectHistory(int steps);

  /**
   * Find possible matrix-element correction object for the dipoles in
   * this state.
   */
  void findMECorrs();

  /** 
   * Returns true if all gluon have an invariant pt above the cutoff.
   */
  bool checkGluonPT();

protected:

  /** @name The virtual functions to be overridden in sub-classes. */
  //@{
  /**
   * Return a simple clone of this object. Should be implemented as
   * <code>return new_ptr(*this);</code> by a derived class.
   */
  virtual ClonePtr clone() const;

  /**
   * Fill the provided set with all pointers to CloneBase objects used
   * in this object.
   */
  inline virtual void fillReferences(CloneSet &) const;

  /**
   * Rebind pointers to other CloneBase objects. Called after a number
   * of interconnected CloneBase objects have been cloned, so that
   * the cloned objects will refer to the cloned copies afterwards.
   *
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   */
  virtual void rebind(const TranslationMap & trans);
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

  /**
   * Determine and return the type of process using information from
   * incoming partons, \a inc, in the sub-process.
   */
  ProcessType getType(const PPair & inc) const;

  /**
   * Determine and store information about the colliding particles
   * using information from incoming partons, \a inc, in the
   * sub-process.
   *
   * @return the pair of momenta of the found incoming particles.
   */
  pair<LorentzMomentum,LorentzMomentum> setIncoming(const PPair & inc);

  /**
   * If the \a partons coming into a hard subprocess goes into an
   * s-channel resonance return this particle, otherwise return null.
   */
  tPPtr singleResonance(const PPair & partons) const;

  /**
   * Given a ColourSinglet object, \a sing, create the corresponding
   * String object together with the necessary Dipole and Parton
   * objects. If \a respectScales is true the created dipoles will get
   * a maximum scale set by the minimum scale of the original partons.
   * If either \a rems is present it is a temporary remant particle
   * which should be handlers separately.
   */
  void createString(const ColourSinglet & sing,
		    bool respectScales = false, tPPair rems = tPPair());

public:

  /** @cond EXCEPTIONCLASSES */
  /**
   * Exception class used if a junction string was found - we don't
   * know how to cascade them yet.
   */
  struct JunctionException: public Exception {};

  /**
   * Exception class used if hadron remnants were found - we can't
   * cascade them yet.
   */
  struct RemnantException: public Exception {};

  /**
   * Exception class used if a string with less than two partons was
   * found.
   */
  struct StringException: public Exception {};
  /** @endcond */

protected:

  /**
   * These are all the emitters included in this state.
   */
  mutable EmitterSet emitters;

  /**
   * These are all the strings included in this state.
   */
  mutable StringSet strings;

private:

  /**
   * A set containing all objects included in this state.
   */
  mutable BaseSet objects;

  /**
   * The Emitter selected to perform the next emission.
   */
  tEmiPtr theSelected;

  /**
   * The pair of remnant partons, if any.
   */
  tRemPair theRemnants;

  /**
   * Information abot the current hard sub-system.
   */
  HardSubSys theHardSubSys;

  /**
   * The process type.
   */
  ProcessType thePType;

  /**
   * The number of emissions made in this state.
   */
  int theNe;

  /**
   * The incoming particles giving the hadronic system.
   */
  tPPair theParticles;

  /**
   * The scattered lepton(s) in case this is a DIS event.
   */
  tPPair theScatteredLeptons;

  /**
   * The scattered quark(s) in case this is a DIS event.
   */
  tPPair theScatteredQuarks;

private:
  /**
   * A Selector containing possible histories to the state.
   * each element will contain a vector with earlier states etc.
   * For each element theSelected is set to the modiefied dipole.
   */
  DipoleStateSelector theHistory;

  /**
   * The selected history.
   */
  DipoleStatePtr theSelectedHistory;

public:

  /**
   * Print out debugging information on std::cerr.
   */
  virtual void debugme() const;

  /**
   * Check integrety of all the emitters. Return false if error is
   * found.
   */
  bool checkIntegrety();

  /**
   * Create and return an index for the given object. Only used in for
   * debugging.
   */
  int index(tcCascadeBasePtr o);

  /**
   * Keep track of indices for debugging purposes.
   */
  ObjectIndexer<int,const Emitter> emiindx;

  /**
   * Keep track of indices for debugging purposes.
   */
  ObjectIndexer<int,const Parton> parindx;

  /**
   * Keep track of indices for debugging purposes.
   */
  ObjectIndexer<int,const String> strindx;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DipoleState> initDipoleState;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleState & operator=(const DipoleState &);

public:

  /** @cond EXCEPTIONCLASSES */
  /** Exception class used by DipoleState if four momentum is
   *  not conserved.
   */
  class MomentumException: public Exception {};

  /** @endcond */

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DipoleState. */
template <>
struct BaseClassTrait<Ariadne::DipoleState,1> {
  /** Typedef of the first base class of DipoleState. */
  typedef Ariadne::CascadeBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DipoleState class and the shared object where it is defined. */
template <>
struct ClassTraits<Ariadne::DipoleState>
  : public ClassTraitsBase<Ariadne::DipoleState> {
  /** Return a platform-independent class name */
  static string className() { return "Ariadne::DipoleState"; }
  /** Return the name of the shared library be loaded to get
   *  access to the DipoleState class and every other class it uses
   *  (except the base class). */
  static string library() { return "libArCascade.so"; }
};

/** @endcond */

}

#include "DipoleState.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DipoleState.tcc"
#endif

#endif /* ARIADNE_DipoleState_H */
