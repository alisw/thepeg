// -*- C++ -*-
#ifndef DIPSY_DiffractiveEventFiller_H
#define DIPSY_DiffractiveEventFiller_H
//
// This is the declaration of the DiffractiveEventFiller class.
//

#include "EventFiller.h"

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

namespace DIPSY {

using namespace ThePEG;

/**
 * The EventFiller class is able to produce an initial
 * ThePEG::Collision from two colliding DipoleStates.
 *
 * @see \ref EventFillerInterfaces "The interfaces"
 * defined for EventFiller.
 */
class DiffractiveEventFiller: public EventFiller {


public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DiffractiveEventFiller();

  /**
   * The destructor.
   */
  virtual ~DiffractiveEventFiller();
  //@}

public:

  /**
   * Exception class for errors in the combinators in the diffractive cascades.
   */
  struct DiffractiveCombinatorics: public Exception {};

  /**
   * Exception class for incorrect cascade structures.
   */
  struct CascadeStructure: public Exception {};

public:

  /** @name Virtual functions which can be overridden in subclasses. */
  //@{
  /**
   * Fill the current collision object in the given DipoleEventHandler
   * with the final state gluons produced when dipole states \a dl and
   * \a dr collides with impact parameters \a b.
   * @return the total interaction probability.
   */

  /**
   * Fill the current collision object in the given DipoleEventHandler
   * with the final state gluons produced when the dipole state \a dl
   * diffractively scatters of the virtual states stored in the eventfiller.
   * @return the total interaction probability.
   */
  virtual double fill(Step & step, DipoleEventHandler & eh, tPPair inc,
			 DipoleState & dl, DipoleState & dr,
			 const ImpactParameters & b) const;

  /**
   * transfer p_minus and p_plus from the elastic state to set the excited state on shell.
   **/
  virtual bool balanceMomenta(DipoleStatePtr excited, DipoleStatePtr elastic) const;

  /**
   * calculates and logs the contribution to the total diffractive cross section.
   **/
  virtual void calculateInclusive(DipoleStatePtr real, CrossSection weight,
				  const ImpactParameters & b) const;

  /**
   * prepares the virtual cascades that all the real states will collide against.
   **/
  virtual void initialiseVirtualCascades(DipoleEventHandler & eh, WaveFunction & WFL,
				  Energy E, double ymax);

  /**
   * randomly selects one of the subcascades of dr with a certain probability for each subcascade
   * and turns dl into that subcascade.
   * returns the weight to be given to the event (0 if no subcascade is chosen)
   **/
  virtual double selectSubcascade(DipoleState & dr, double selectionProbability) const;

  /**
   * randomly selects one of the partons that is not on-shell
   * and puts that parton and its parents on shell in the state.
   **/
  virtual double addRandomChain(DipoleState & dr) const;

  /**
   * randomly selects one of the single-interaction subcascades of dr with
   * a certain probability for each subcascade
   * and turns dl into that subcascade.
   * returns the weight to be given to the event (0 if no subcascade is chosen)
   **/
  virtual double selectSingleSubcascade(DipoleState & dr) const;

  /**
   * randomly selects one of the subcascades of dr with N interactions
   * and turns dl into that subcascade.
   * returns an approximation of the number of possible choices
   **/
  virtual double selectNSubcascade(DipoleState & dr, int N) const;

  /**
   * coutns the number of subcacscades with exactly one interaction
   **/
  virtual int countSingleStates(DipoleState & dr) const;

  /**
   * counts and returns the number of subcascades of the state
   **/
  virtual int nSubcascades(DipoleState & dr) const;

  /**
   * changes the state to its substate of choice.
   **/
  virtual double pickSubcascade(DipoleState & dr, int choice) const;

  /**
   * changes the state to its substate o.
   **/
  virtual void pickSubcascadeFromInteraction(DipoleState & dr, PartonPtr p) const;

  /**
   * adds a parton and its parents to the on-shell state.
   **/
  // virtual void addChainFromInteraction(DipoleState & dr, PartonPtr p) const;

  /**
   * adds the two partons of the dipole to the set x.
   **/
  virtual void addPartonsToSet(DipolePtr d, set<PartonPtr> x) const;

  /**
   * Set the parton and all its ancestors on shell.
   **/
  virtual void recursiveSetOnShell(PartonPtr p) const;

  /**
   * return more likely cascades and removes(return 0.0) less likely ones.
   * A state with small dipole will probably get removed, but the ones that doesnt are returned with
   * a larger weight.
   **/
  double reweightCascade(DipoleState & dr, vector<PartonPtr> chainCaps) const;

  /**
   * return more likely cascades and removes(return 0.0) less likely ones.
   * A state with small dipole will probably get removed, but the ones that doesnt are returned with
   * a larger weight.
   **/
  double estimateF(DipoleStatePtr dr) const;

  /**
   * sets all the partons in the state to off-shell
   **/
  virtual void setOffShell(DipoleState & dr) const;

  /**
   * sets all the valence partons in the state on-shell.
   **/
  virtual void setValenceOnShell(DipoleState & dr) const;

  /**
   * sets the two partons in the dipole on shell.
   **/
  virtual void setPartonsOnShell(DipolePtr dip) const;

  /**
   * changes the subcascade from the dipole to its substate of choice.
   **/
  virtual double pickSubcascade(DipolePtr dip, int choice) const;

  /**
   * checks and fixes the parent structure so that the on shell partons only
   * are connected to each other. The others are bypassed.
   * This may leave some partons without children if the swing is present.
   **/
  virtual void updateParentStructure(DipoleStatePtr dr) const;

  /**
   * finds the nonvalence partons that does not have any children.
   * These are at the end of the chains, and can loosely be interpreted as the interactions.
   **/
  virtual vector<PartonPtr> childlessChildren(DipoleStatePtr dr) const;

  /**
   * removes all valence partons from the provided vector
   **/
  virtual void removeValence(vector<PartonPtr> partons) const;

  /**
   * counts and returns the number of subcascades from this dipole
   **/
  virtual int nSubcascades(DipolePtr dip) const;

  /**
   * makes dl real and on-shell and does all the reweighting and ordering bussiness.
   * p- (and pt) is balances by dr, but dr is otherwise not touched.
   **/
  virtual void makeReal(DipoleState & dl, vector<PartonPtr> chainCaps, DipoleState & dr) const;

  /**
   * Calculates the amplitude for the real state dr to diffractively scatter on
   * the virtual cascades stored in the eventfiller.
   **/
  virtual double amplitude(DipoleState & dr, vector<PartonPtr> chainCaps,
					 const ImpactParameters & b) const;

  /**
   * Calculates the elastic amplitude for dr on the virtual cascades stored in the eventfiller.
   **/
  virtual double elasticAmplitude(DipoleState & dr, const ImpactParameters & b) const;

  /**
   * removes the chain caps from the state dr.
   * Returns true if there is still at least one dipole left.
   **/
  virtual DipoleStatePtr steriliseCaps(DipoleStatePtr dr) const;

  /**
   * returns a DipoleState with the two neighbouring dipoles to each parton in the vector.
   **/
  virtual vector<DipoleStatePtr> extractI(vector<PartonPtr> chainCaps) const;

  /**
   * returns a DipoleState with the dipole that emitted each parton in the vector.
   **/
  virtual vector<DipoleStatePtr> extractH(vector<PartonPtr> chainCaps) const;

  /**
   * returns the original emission rapidity for each parton in the vector.
   **/
  virtual vector<double> extractIntY(vector<PartonPtr> chainCaps) const;

  /**
   * returns a string with a single proton with the 3-momentum of the original state
   **/
  virtual DipoleState::String makeProton(DipoleStatePtr) const;


  //inline functions

  /**
   * get the mode.
   */
  inline int nElasticCascades() const {
    return theNElasticCascades;
  }

  /**
   * set the mode.
   */
  inline void nElasticCascades(int x) {
    theNElasticCascades = x;
  }

  /**
   * get the minimum rapidity gap.
   */
  inline double minRapGap() const {
    return theMinRapGap;
  }

  /**
   * set the minimum rapidity gap.
   */
  inline void minRapGap(double x) {
    theMinRapGap = x;
  }

  /**
   * get the number of interactions.
   */
  inline int nInteractions() const {
    return theNInteractions;
  }

  /**
   * set the number of interactions.
   */
  inline void nInteractions(int x) {
    theNInteractions = x;
  }

  /**
   * get the maximum rapidity gap.
   */
  inline double maxRapGap() const {
    return theMaxRapGap;
  }

  /**
   * set the maximum rapidity gap.
   */
  inline void maxRapGap(double x) {
    theMaxRapGap = x;
  }

  /**
   * get the minimum rapidity.
   */
  inline double minY() const {
    return theMinY;
  }

  /**
   * set the minimum rapidity.
   */
  inline void minY(double x) {
    theMinY = x;
  }

  /**
   * get the maximum rapidity.
   */
  inline double maxY() const {
    return theMaxY;
  }

  /**
   * set the maximum rapidity.
   */
  inline void maxY(double x) {
    theMaxY = x;
  }

  /**
   * reset the number of subcascades.
   */
  inline void resetNSubCascades() {
    theNSubCascades = 0.0;
  }

  /**
   * increase the number of subcascades.
   */
  inline void increaseNSubCascades(double x) const {
    theNSubCascades += x;
  }

  /**
   * reset the number of subcascades.
   */
  inline double nSubCascades() {
    return theNSubCascades;
  }

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
   * Is this where this should be declared?? :o
   */
  virtual void dofinish();

private:
  mutable CrossSection totalXSec;
  mutable QTY<4,0,0>::Type totalSqr;
  mutable int neve;

private:

  //The linera combination of virtual cascades from the elastic particle
  vector<DipoleStatePtr> virtualCascades;

  //the number of pregenerated elstic cascades
  int theNElasticCascades;

  //the minimum rapidity gap
  double theMinRapGap;

  //the number of interactions
  int theNInteractions;

  //the total number of subcascades encountered
  mutable double theNSubCascades;

  //the maximum rapidity gap
  double theMaxRapGap;

  //the minimum rapidity that every excited state cross.
  double theMinY;

  //the maximum rapidity that no excited state cross.
  double theMaxY;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DiffractiveEventFiller & operator=(const DiffractiveEventFiller &);

};

}

#endif /* DIPSY_DiffractiveEventFiller_H */
