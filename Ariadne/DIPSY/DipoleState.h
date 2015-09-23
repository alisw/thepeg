// -*- C++ -*-
#ifndef DIPSY_DipoleState_H
#define DIPSY_DipoleState_H
//
// This is the declaration of the DipoleState class.
//

#include "ThePEG/Config/ThePEG.h"
#include "DipoleState.fh"
#include "Dipole.h"
#include "DipoleEventHandler.fh"
#include "WFInfo.h"
#include "DipoleXSec.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ShadowParton.h"

#include "ThePEG/Analysis/FactoryBase.h"

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the DipoleState class.
 */
class DipoleState: public PersistentBase {

public:

  /**
   * Copy FList typedef from DipoleXSec.
   */
  typedef DipoleXSec::FList FList;

  /**
   * A String is simply a vector of colour-connected partons.
   */
  typedef vector<PartonPtr> String;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The standard constructor taking the controlling
   * DipoleEventHandler as argument.
   */
  inline DipoleState(const DipoleEventHandler & h, WFInfoPtr wf = WFInfoPtr())
    : thePlus(ZERO), theMinus(ZERO), theMinusDeficit(ZERO), theHandler(&h),
      theWFInfo(wf), theWeight(1.0), doTakeHistory(false), theYmax(0.0),
      theCollidingEnergy(ZERO) {}

  /**
   * The default constructor.
   */
  inline DipoleState()
    : thePlus(ZERO), theMinus(ZERO), theMinusDeficit(ZERO), theWeight(1.0),
      doTakeHistory(false), theYmax(0.0), theCollidingEnergy(ZERO)  {}

  /**
   * The copy constructor.
   */
  DipoleState(const DipoleState &);

  /**
   * The destructor.
   */
  virtual ~DipoleState();
  //@}

public:

  /**
   * Runs through a final state evolution checking for swings only.
   */
  void swingFS(double ymin, double ymax);

  /**
   * Does some colour recconection for the valence partons to account for
   * original 3 colour colour structure of a proton.
   */
  void normaliseValenceCharge(int mode);

  /**
   * Makes sure the non-participating nucleons of an AA collision have the
   * colour flow of the original triangles.
   */
  void restoreNonparticipants();

  /**
   * counts and returns the number of spectating nucleons.
   */
  int numberOfSpectators() const;

  /**
   * Makes sure the the two partons are conected by a dipole with p1 as first parton.
   * returns false if it fails, and they are left not connected.
   */
  bool restoreDipole(PartonPtr p1, PartonPtr p2);

  /**
   * Evolve this state from the given minimum to the given maximum
   * rapidity, assuming there will be p- coming from the other state.
   */
  void evolve(double ymin, double ymax);

  /**
   * Makes the state merge with another colliding state.
   */
  DipoleStatePtr collide(DipoleStatePtr otherState,
			 const vector<FList::const_iterator> & sel,
                         const ImpactParameters & b);

  /**
   * Makes the state merge with another colliding state.
   */
  DipoleStatePtr merge(DipoleStatePtr otherState);

  /**
   * mirrors the state in rapidity around y0.
   */
  void mirror(double y0);

  /**
   * moves all partons and their pT according to the imparctparameter b.
   */
  void translate(const ImpactParameters & b );

  /**
   * Fix up valens partons if they were not of the correct flavour or
   * if they should collapse into a hadron.
   */
  virtual void fixValence(Step & step) const;

  /**
   * returns the positions of all partons in the state,
   * together with the size of the largest neighbouring dipole.
   */
  vector<pair<Parton::Point, InvEnergy> > points();

  /**
   * Checks through all partons and reabsorbs the ones which are
   * considered enough gain in some kind of virtuality.
   */
  void absorbSmallDipoles();

 /**
   * Removes all the parents of the dipoles, and sets the currently
   * active dipoles as the initial ones.
   **/
  void makeOriginal();

  /**
   * returns teh average distance in rapidity of the interacting dipoles.
   */
  double avYInInt() const;

  /**
   * Controls and couts various info about the state. Returns false
   * if something is too wrong.
   */
  bool diagnosis(bool print) const;

  /**
   * Check energy-momentum conservation.
   */
  void checkFSMomentum() const;

  /**
   * Check energy-momentum conservation.
   */
  void checkFSMomentum(const Step & step) const;

  /**
   * Prints the state to file.
   */
  void saveGluonsToFile(double weight) const;

  /**
   * Returns the partons sorted in Strings.
   */
  vector<DipoleState::String> strings();

  /**
   * Returns the active partons.
   */
  list<PartonPtr> getPartons() const;

  /**
   * Returns the active dipoles.
   */
  list<DipolePtr> getDipoles() const;

  /**
   * Returns all dipoles.
   */
  set<DipolePtr> getAllDipoles() const {
    return allDipoles;
  }

  /**
   * Makes the state reabsorb a parton.
   */
  void reabsorb(PartonPtr p);

  /**
   * Makes the state reabsorb a parton, and continues to absorb as
   * long as the dipoles gets smaller along the colour chain.
   */
  void recursiveReabsorb(PartonPtr p);

  /**
   * Forces one of the dipoles connected to the parton 
   * to immediately swing. But not with higher rapdity
   * step than /a ymax1 of same colour, or /a ymax2 of
   * another colour if same colour cant be found.
   * /a ymax = 0 mean no limit.
   */
  bool forceSwing(PartonPtr d, double ymax1, double ymax2);

  /**
   * Shuffles p+ and p- between the left and right part of the state
   * so that the particles are really on shell. Assumes m=0 atm.
   */
  void balanceMomenta();

  /**
   * Check if shadows has been setup.
   */
  bool hasShadows() const {
    return !theIncomingShadows.empty();
  }

  /**
   * Set up shadow partons for all valence partons, thus initializing
   * the shadow parton procedure.
   */
  void setupShadows();

  /**
   * Reset all interaction flags and interactions in all shadows.
   */
  void resetShadows();

  /**
   * Reset all interaction flags in all shadows.
   */
  void resetInteractedShadows();

  /**
   * Return the incoming momentum of the given valens shadow parton.
   */
  LorentzMomentum incomingMomentum(tcSPartonPtr valence, int mode);

  /**
   * Get the weight associated with the generation of this dipole state.
   */
  inline double weight() const {
    return theWeight;
  }

  /**
   * Set the weight associated with the generation of this dipole state.
   */
  inline void weight(double x) {
    theWeight = x;
  }

  /**
   * Get the p+ that the original particle brought.
   */
  inline Energy plus() const {
    return thePlus;
  }

  /**
   * Set the p+ that the original particle brought.
   */
  inline void plus(Energy E) {
    thePlus = E;
  }

  /**
   * Get the p- that the original particle brought.
   */
  inline Energy minus() const {
    return theMinus;
  }

  /**
   * Set the p- that the original particle brought.
   */
  inline void minus(Energy E) {
    theMinus = E;
  }

  /**
   * Get the p- that the original particle was missing.
   */
  inline Energy minusDeficit() const {
    return theMinusDeficit;
  }

  /**
   * Set if the state should save its history or not.
   */
  inline void takeHistory(bool b) {
    doTakeHistory = b;
  }

  /**
   * saves the current state into the states history, but only if
   * doTakeHistory is true.
   */
  inline void save() {
    if ( doTakeHistory ) theHistory.push_back( clone() );
  }

  /**
   * The list of initial dipoles.
   */
  inline const vector<DipolePtr> & initialDipoles() const {
    return theInitialDipoles;
  }

  /**
   * Adds a dipole to the initial Dipoles.
   * Mainly for testing purposes.
   */
  inline void addDipole(Dipole & dip) {
    theInitialDipoles.push_back(& dip);
  }

  /**
   * Get additional info about the wavefunction used to create this state.
   */
  inline WFInfoPtr WFInfo() const {
    return theWFInfo;
  }

  /**
   * Get additional info about the wavefunction used to create this state.
   */
  inline const WaveFunction & wf() const {
    return WFInfo()->wf();
  }

  /**
   * Set additional info about the wavefunction used to create this state.
   */
  inline void WFInfo(WFInfoPtr x) {
    theWFInfo = x;
  }

  /**
   * The controlling DipoleEventHandler.
   */
  inline const DipoleEventHandler & handler() const {
    return *theHandler;
  }

  /**
   * Set the controlling DipoleEventHandler.
   * added by CF to access from emitter.
   */
  inline void handler(tcDipoleEventHandlerPtr hdl) {
    theHandler = hdl;
  }

  /**
   * Create a Dipole belonging to this state.
   */
  inline tDipolePtr createDipole() {
    DipolePtr d = new_ptr(Dipole());
    d->theDipoleState = this;
    allDipoles.insert(d);
    return d;
  }

  /**
   * Generate a consistent colour index for the given Dipole.
   */
  void generateColourIndex(tDipolePtr);

  /**
   * Make all colours belong to the same system.
   */
  void unifyColourSystems(int isys = 0);
  /**
   * Returns all active dipoles of a certain colour. Optionally only
   * return dipoles which have been touched since the last emission.
   */
  inline const vector<tDipolePtr> &
  swingCandidates(int c, bool touched = false) const {
    static vector<tDipolePtr> dummy;
    if ( !touched ) return theSwingCandidates[c];
    if ( c >= int(theTouchedSwingCandidates.size()) ) return dummy;
    return theTouchedSwingCandidates[c];
  }

  /**
   * adds the dipole, and all its children, to theSwingCandidates list.
   * includes interacting dipoles.
   */
  void sortDipoleFS(Dipole &);

  /**
   * Sort in colour all the dipoles originating from initialDipoles.
   * includes interacting dipoles.
   */
  void sortDipolesFS();

  /**
   * Sort in colour all the dipoles originating from initialDipoles.
   * includes interacting dipoles.
   */
  void sortFSDipoles();

  /**
   * Return the total lambda measure for all final state dipoles.
   */
  pair<double,int> lambdaMeasure(Energy2 scale = 1.0*GeV2,
				 FactoryBase::tH1DPtr histlength = 0,
				 FactoryBase::tH1DPtr histmass = 0) const;

  /**
   * adds the dipole, and all its children, to theSwingCandidates list.
   */
  void sortDipole(Dipole &);

  /**
   * Sort in colour all the dipoles originating from initialDipoles.
   */
  void sortDipoles();

  /**
   * Returns a list of all the active nonvalence partons that has not interacted.
   */
  const list<PartonPtr> virtualPartons() const;

  /**
   * Calculates and returns the set of the loops in the state.
   */
  const set< list<PartonPtr> > loops() const;

  /**
   * Returns the highest rapidity in any parton in the state.
   */
  const double highestY() const;

  /**
   * Returns the highest rapidity in any parton in the state.
   */
  const double lowestY() const;

  /**
   * Returns the rapidity the state has been evolved to.
   */
  const double ymax() const;

  /**
   * Gets the energy the colliding state will supply.
   */
  const Energy collidingEnergy() const;

  /**
   * Sets the energy the colliding state will supply.
   */
  void collidingEnergy(Energy);

  /**
   * Return the incoming particle(s) mapped to the corresponding valence particles.
   */
  inline const map<PPtr, vector<PartonPtr> > & incoming() const {
    return theIncoming;
  }

  /**
   * Return a ThePEG::Particle produced from the given
   * DIPSY::Parton. I nocreate is true, only return a particle which
   * has previously been created, do not create any new ones.
   */
  tPPtr getParticle(tcPartonPtr parton, bool nocreate = false) const;

  /**
   * Given an output iterator fill the corresponding container with
   * final, undecayed dipoles. This is done recursively, ie. if this
   * dipole has decayed, the extract methods of the two children are
   * called instead.
   */
  template <typename OutputIterator>
  inline void extract(OutputIterator it) const {
    for ( int i = 0, N = initialDipoles().size(); i < N; ++i )
      initialDipoles()[i]->extract(it);
  }

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
   * Clone this state and all the dipoles and partons in it.
   */
  DipoleStatePtr clone();

protected:

  /**
   * The p+ and p- that the original particle had.
   */
  Energy thePlus;
  Energy theMinus;

  /**
   * Keeps track of extra missing p- that is not reflected in the
   * valence partons.
   */
  Energy theMinusDeficit;

  /**
   * The controlling DipoleEventHandler.
   */
  tcDipoleEventHandlerPtr theHandler;

  /**
   * The list of initial dipoles.
   */
  vector<DipolePtr> theInitialDipoles;

  /**
   * The active dipoles sorted by colour. /CF
   */
  vector< vector<tDipolePtr> > theSwingCandidates;

  /**
   * The active dipoles which have been touched since the last
   * emission, sorted by colour
   */
  vector< vector<tDipolePtr> > theTouchedSwingCandidates;

  /**
   * Additional info about the wavefunction used to create this state.
   */
  WFInfoPtr theWFInfo;

  /**
   * The weight associated with the generation of this dipole state.
   */
  double theWeight;

  /**
   * If the history should be saved and displayed or not.
   */
  bool doTakeHistory;

  /**
   * The rapidity the state has been evolved to.
   */
  double theYmax;

  /**
   * How much energy the meeting particle have.
   */
  Energy theCollidingEnergy;

  /**
   * The history of this state.
   */
  vector<DipoleStatePtr> theHistory;

  /**
   * The set of all dipoles belonging to this state.
   */
  set<DipolePtr> allDipoles;

  /**
   * A map relating incoming particles with the valens partons.
   */
  map<PPtr, vector<PartonPtr> > theIncoming;

  /**
   * A map relating a produced ThePEG::Particle to its DIPSY::Parton.
   */
  mutable map<tcPartonPtr,PPtr> theProducedParticles;

  /**
   * The shadowing incoming partons.
   */
  vector<SPartonPtr> theIncomingShadows;

  /**
   * The propagator momenta of each valence shadow. Set once and for
   * all as the total momentum of the state minus all other valence
   * shadows.
   */
  map<tcSPartonPtr,LorentzMomentum> theShadowPropagators;

protected:

  /**
   * Exception class for badly connected dipoles.
   */
  struct DipoleConnectionException: public Exception {};

  /**
   * Exception class for bad kinematics.
   */
  struct DipoleKinematicsException: public Exception {};

  /**
   * Exception class for bad DGLAP checks.
   */
  struct DipoleDGLAPSafeException: public Exception {};

public:

  /**
   * Helper function for fixValence implementation. Take the momentum
   * \a p and return a new momentum with the mass \a m, such that the
   * total rest massof the system when combined with a fererence
   * momentum, \a ref, and the reference momentum itself is
   * preserved. Note that this will change the total momentum of the
   * system, but this is fixed in the end by a simple global boost in
   * the fixValence method. If \a Rshift is non-null, the
   * corresponding object is set to the boost from the old total frame
   * to the current one.
   */
  static LorentzMomentum changeMass(LorentzMomentum p, Energy m, LorentzMomentum ref,
				    LorentzRotation * Rshift = 0);

  /**
   * Printout the shadow structure if available.
   */
  void debugShadowTree() const;

  void checkShadowMomentum(Sum20Momentum & sum20,
			   const ImpactParameters * b = 0) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleState & operator=(const DipoleState &);

};

}

#endif /* DIPSY_DipoleState_H */
