// -*- C++ -*-
#ifndef Ariadne5_DipoleState_H
#define Ariadne5_DipoleState_H
//
// This is the declaration of the DipoleState class.
//

#include "Ariadne/Cascade/CascadeBase.h"
#include "DipoleState.fh"
#include "DipoleBase.h"
#include "Parton.h"
#include "EmissionGenerator.h"
#include "Emission.h"
#include "RemnantParton.h"
#include "Resonance.fh"
#include "ThePEG/Utilities/ObjectIndexer.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The DipoleState class describes a complete state of partons and
 * dipoles in a Ariadne cascade at a given stage of the process.
 */
class DipoleState: public CascadeBase {

public:

  /**
   * A vector of CascadeBase pointers.
   */
  typedef set<CascadeBasePtr> BaseSet;

  /**
   * A set of dipole pointers.
   */
  typedef list<tDBPtr> DipoleSet;

  /**
   * The SaveDipoleState helper class is a friend.
   */
  friend class SaveDipoleState;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor. Optional argument is the incoming
   * particles to the collision.
   */
  DipoleState(tPPair inc = tPPair());

  /**
   * The destructor.
   */
  virtual ~DipoleState();
  //@}

public:

  /**
   * Setup this DipoleState from a SubProcess.
   */
  void setup(SubProcess &);

  /**
   * Setup this DipoleState from a simple set of final-state particles
   */
  void setup(const set<tPPtr> &);

  /**
   * Select a dipole. Go through the list of active DipoleBase objects
   * and tell them to calculate the scale of the next emission and
   * select the Emitter with largest transverse momentum. Subsequent
   * calls to selected() will return the selected dipole.
   *
   * @param rhomin the minimum allowed evolution scale.
   *
   * @param pt2max the maximum allowed evolution scale.
   *
   * @return the generated evolution scale. If less that rhomin, no
   * emission was generated and subsequent calls to selected() will
   * return null.
   */
  Energy select(Energy rhomin, Energy rhomax);

  /**
   * Perform an emission.
   *
   * @return false if no emission was possible.
   */
  bool perform();

  /**
   * After the cascade is done, fill the produced partons in the
   * supplied Step.
   */
  void fill(Step &) const;

  /** @name Simple access functions. */
  //@{
  /**
   * The next Emission
   */
  inline tEmPtr selected() const {
    return theSelected;
  }

  /**
   * Set the next Emission.
   */
  inline void selected(const tEmPtr & x) {
    theSelected = x;
  }

  /**
   * The number of emissions made from this state so far.
   */
  inline int nEmissions() const {
    return theNe;
  }

  /**
   * Return the colliding particles if present.
   */
  tPPair incoming() const {
    return theIncoming;
  }

  /**
   * Return the final-state partons. This includes all final-state
   * particles, also remnants.
   */
  const set<tParPtr> & finalState() const {
    return theFinalState;
  }

  /**
   * Add a parton to the final state.
   */
  void addFS(tParPtr p) {
    theFinalState.insert(p);
  }

  /**
   * Add a parton to the hard final state.
   */
  void addHardFS(tParPtr p) {
    theHardFinalState.insert(p);
    addFS(p);
  }

  /**
   * Add a parton to the hadronic final state.
   */
  void addHadronicFS(tParPtr p) {
    theHadronicFinalState.insert(p);
    addHardFS(p);
  }

  /**
   * Remove the given parton from the final state and forget it ever
   * existed.
   */
  void forgetParton(tParPtr p);

  /**
   * Remove the given dipole from the active dipoles and forget it ever
   * existed.
   */
  void forgetDipole(tDBPtr d);

  /**
   * Return the partons in the hard final state. Returns finalState()
   * excluding all soft remnants.
   */
  const set<tParPtr> & hardFS() const {
    return theHardFinalState;
  }

  /**
   * Return the partons in the hadronic final state excluding all soft
   * and non-coloured remnants.
   */
  const set<tParPtr> & hadronicFS() const {
    return theHadronicFinalState;
  }

  /**
   * Return the remnants. Are a subset of theFinalState
   */
  const pair< vector<tRemParPtr>, vector<tRemParPtr> > & remnants() const {
    return theRemnants;
  }

  /**
   * Return the set of active dipoles.
   */
  const set<tDBPtr> & activeDipoles() const {
    return active;
  }

  /**
   * Return massive resonances.
   */
  const vector<tResPtr> & resonances() const {
    return theResonances;
  }
  //@}

public:

  /** @name Functions relating to the book-keeping of included objects. */
  //@{
  /**
   * Create a new object to be included in this DipoleState.
   */
  template <typename Class>
  inline typename Ptr<Class>::pointer create() {
    typedef typename Ptr<Class>::pointer PTR;
    PTR obj = ptr_new<PTR>();
    objects.insert(objects.end(), obj);
    obj->state(this);
    if ( tDBPtr d = dynamic_ptr_cast<tDBPtr>(obj) ) {
      active.insert(d);
      dipindx(d);
    }
    if ( tParPtr p = dynamic_ptr_cast<tParPtr>(obj) ) parindx(p);
    return obj;
  }

  /**
   * Create a Parton of the given type \a pd. Optionally inherit the
   * system properties of the \a parent parton and insert in the
   * hadronic final state if \a hfs is true.
   */
  tParPtr create(tcPDPtr pd, tcParPtr parent = tcParPtr(), bool hfs = true);

  /**
   * Remove a Parton from the final state. (The parton will still
   * exist in the DipoleState).
   */
  void remove(tParPtr p) {
    theFinalState.erase(p);
    theHardFinalState.erase(p);
    theHadronicFinalState.erase(p);    
  }

  /**
   * Clone this DipoleState. Also cloning all included objects and fixing up
   * their inter-dependence.
   */
  DipoleStatePtr fullclone() const;

  /**
   * Clone this DipoleState. Also clone all included objects but do
   * not fix up their inter-dependence, ie. the cloned objects in the
   * new DipoleState may still point to objects in the this. Before
   * the new DipoleState can be used, its postclone() function must
   * be called with the \a trans object as argument. The \a trans
   * object will contain the translation map between the objects in
   * the old and new DipoleState.
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
   * Calculate the total momentum of the dipole state. 
   */
  const Lorentz5Momentum & sumTotalMomentum();

  /**
   * Return the total momentum of the dipole state. This is
   * recalculated after each emission.
   */
  const Lorentz5Momentum & totalMomentum() const {
    return theTotalMomentum;
  }

  /**
   * Calculate the total momentum of the hardFinalState(). 
   */
  const Lorentz5Momentum & sumHardMomentum();

  /**
   * Return the total momentum of the hardFS(). This is
   * recalculated after each emission.
   */
  const Lorentz5Momentum & hardMomentum() const {
    return theHardMomentum;
  }

  /**
   * Return the total momentum of the hadronicFS(). This is
   * recalculated after each emission.
   */
  const Lorentz5Momentum & hadronicMomentum() const {
    return theHadronicMomentum;
  }

  /**
   * Return the total sting length in terms of the lambda measure
   * using the given scale.
   */
  pair<double,int> lambdaMeasure(Energy2 scale = 1.0*GeV2) const;

  /**
   * Return thenumber of dipoles which stretches across the 0-rapidity line.
   */
  int crossings() const;

  /**
   * Return the number of gluons below cutoff (for debugging).
   */
  int gluonsBelow() const;

  /**
   * Return the length in true rapidity of all stings divided by the
   * total rapidity span of the state (for debugging).
   */
  double folding() const;

  /**
   * Transform the hadronicFS().
   */
  void transformHadronicState(const LorentzRotation & R);

  /**
   * Touch all remnants to indicate that the hadronicFS() has changed
   * collectively.
   */
  void touchHadronicState();

  /**
   * Untouch all remnants to indicate that the hadronicFS() has
   * reverted collectively.
   */
  void untouchHadronicState();

  /**
   * Print out debugging information on std::cerr.
   */
  virtual void debugme() const;

  /**
   * Print out more debugging information on std::cerr. 
   */
  void debugEmissions() const;

  /**
   * Check integrity of all dipoles. Return false if error is
   * found.
   */
  bool checkIntegrity();

  /**
   * Create and return an index for the given object. Only used in for
   * debugging.
   */
  int index(tcCascadeBasePtr o) const;

  /**
   * For the given particle set, remove all particles that have
   * decayed and replace them with their children.
   */
  static void fixupDecays(set<tPPtr> & out);

  /**
   * Purge all gluons with transvese momentum less that the given cut.
   */
  void purgeGluons(Energy cut);

  /**
   * Purge the outgoing gluon in the given dipole.
   */
  void purgeGluon(QCDDipole & d);

  /**
   * For the given particle insert recursively all undecayed decay
   * products in the given vector.
   */
  static bool addChildren(tPVector & v, tPPtr p);

public:

  /**
   * Trace a dipole to find (a piece of) string.
   */
  static pair<tcQCDPtr,tcQCDPtr> StringEnds(tcQCDPtr);

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
  virtual void fillReferences(CloneSet &) const;

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

private:

  /**
   * The SubProcess from which this state was created (may be null).
   */
  tSubProPtr subprocess;

  /**
   * The colliding particles in the SubProcess or 0 if none was
   * present.
   */
  tPPair theIncoming;

  /**
   * The set of all diferent objects belonging to this state.
   */
  BaseSet objects;

  /**
   * The final state partons.
   */
  set<tParPtr> theFinalState;

  /**
   * The hard final state partons, excluding all soft remnants.
   */
  set<tParPtr> theHardFinalState;

  /**
   * The hadronic final state, excluding all soft and non-coloured
   * remnants.
   */
  set<tParPtr> theHadronicFinalState;

  /**
   * The remnants. Are a subset of theFinalState
   */
  pair< vector<tRemParPtr>, vector<tRemParPtr> > theRemnants;

  /**
   * The massive resonances.
   */
  vector<tResPtr> theResonances;

  /**
   * The total momentum of this state.
   */
  Lorentz5Momentum theTotalMomentum;

  /**
   * The total momentum of the hardFS().
   */
  Lorentz5Momentum theHardMomentum;

  /**
   * The total momentum of the hadronicFS().
   */
  Lorentz5Momentum theHadronicMomentum;

  /**
   * The active Dipoles.
   */
  set<tDBPtr> active;

  /**
   * The set of EmissionGenerators corresponding to the dipoles
   * included in this state.
   */
  set<EmissionGenerator> generators;

  /**
   * The number of emissions made in this state.
   */
  int theNe;

  /**
   * The Emitter selected to perform the next emission.
   */
  EmPtr theSelected;

  /**
   * The list of emissions treated in this DipoleState.
   */
  vector<EmPtr> emissions;

  /**
   * Keep track of indices for debugging purposes.
   */
  mutable ObjectIndexer<int,const DipoleBase> dipindx;

  /**
   * Keep track of indices for debugging purposes.
   */
  mutable ObjectIndexer<int,const Parton> parindx;

private:

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

  /** Exception class used by DipoleState if the SubProcess was
   *  querried when none was available.
   */
  class SubProcessException: public Exception {};
  /** @endcond */

};

/**
 * Helper class to be able to revert changes made to a DipoleState.
 */
class SaveDipoleState {

public:

  /**
   * The constructor taking a pointer to the \a state to be saved. If
   * necessary (or if \a force is true) the state will be cloned.
   */
  SaveDipoleState(tDipoleStatePtr state, bool force = false);

  /**
   * Return the translation map relating the objects in the original
   * state with the ones in the backup state.
   */
  const DipoleState::TranslationMap & translationMap() const {
    return trans;
  }

  /**
   * Return a copy of the original state.
   */
  tDipoleStatePtr revert();

private:

  /**
   * The backup state.
   */
  DipoleStatePtr backup;

  /**
   * The translation map relating the objects in the original state
   * with the ones in the backup state.
   */
  DipoleState::TranslationMap trans;

  /**
   * Flag to indicate that precloning was forced.
   */
  bool forced;

};

}

#endif /* Ariadne5_DipoleState_H */
