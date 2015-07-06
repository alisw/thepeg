// -*- C++ -*-
#ifndef DIPSY_RealParton_H
#define DIPSY_RealParton_H
//
// This is the declaration of the RealParton class.
//

#include "ThePEG/Config/ThePEG.h"
#include "RealParton.fh"
#include "RealPartonState.fh"
#include "Parton.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the RealParton class.
 */
class RealParton: public Base {

public:

  /**
   * Internal class to define the rapidity ordering of RealParton
   * objects.
   */
  struct RealPartonOrder {
    /** The actual ordering function. */
    bool operator()(tRealPartonPtr a, tRealPartonPtr b) const {
      return a->theParton->oY() < b->theParton->oY();
    }
  };

  /**
   * A rapidity-ordered set of RealParton objects.
   */
  typedef set<tRealPartonPtr,RealPartonOrder> RealPartonSet;

  /**
   * A recoil with a parton, and how much momenta was TAKEN from that parton.
   */
  typedef pair<tRealPartonPtr, pair<TransverseMomentum, Energy> > Recoil;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  RealParton(tPartonPtr);

  /**
   * The destructor.
   */
  virtual ~RealParton();
  //@}

public:

  /**
   * The real state it belongs to.
   */
  tRealPartonStatePtr realState;

  /**
   * The possible verdicts for a partons existance.
   */
  enum status {YES, NO};

  /**
   * If this parton should be part of the final state.
   */
  status keep;
  status cKeep;

  /**
   * If this parton is interacting, and with what dipole.
   */
  tDipolePtr interacting;
  tDipolePtr cInteracting;

  /**
   * If the parton has an interaction on the first and second side in colour flow.
   **/
  bool firstInt;
  bool secondInt;
  bool cFirstInt;
  bool cSecondInt;

  /**
   * The number of mothers the parton has.
   **/
  int nMothers;

  /**
   * The mothers in the final state.
   */
  pair<tRealPartonPtr,tRealPartonPtr> mothers;

  /**
   * The original mothers from the evolution.
   */
  pair<tRealPartonPtr,tRealPartonPtr> oMothers;

  /**
   * The last pair of mothers that has been part of a consistent state.
   */
  pair<tRealPartonPtr,tRealPartonPtr> cMothers;

  /**
   * The mother in the final state, in the model with a single mother.
   */
  tRealPartonPtr mother;
  tRealPartonPtr cMother;
  tRealPartonPtr oMother;

  /**
   * The children emitted on the first mothers and second mothers side
   * respectively.
   */
  pair<RealPartonSet, RealPartonSet> children;
  pair<RealPartonSet, RealPartonSet> oChildren;
  pair<RealPartonSet, RealPartonSet> cChildren;

  /**
   * The parton.
   */
  tPartonPtr theParton;

  /**
   * The interactions this parton is connected to.
   */
  set<tDipolePtr> interactions;
  set<tDipolePtr> cInteractions;

  /**
   * The distance to the parton it will be connected to in an interaction,
   * if the parton is interacting.
   */
  InvEnergy intDist;
  InvEnergy cIntDist;

  /**
   * Dipoles towards mothers that have been checked for DGLAP supression,
   * and the scale needed to resolve the dipole.
   */
  map<tRealPartonPtr, InvEnergy2> resolutionScales;

  /**
   * The momentum of the particle.
   */
  TransverseMomentum pT;
  Energy plus;
  Energy minus;
  double y;

  /**
   * The plus the parton took from the first and second effective mother respectively.
   */
  Energy plus1, plus2;

  /**
   * The transverse momenta and p+ that has been taken from partons at emission.
   **/
  list< Recoil > recoils;

  /**
   * The transverse momentum and p+ gotten from the colliding state.
   **/
  pair< TransverseMomentum, Energy > interactionRecoil;
  pair< TransverseMomentum, Energy > cInteractionRecoil;

  /**
   * The future emissions that have effected this parton.
   **/
  RealPartonSet future;

  /**
   * the set of partons that will merge into one after interaction.
   * Points to a set in the realpartonstate.
   **/
  int fluct;
  int cFluct;

  /**
   * The recoilSwings that were made when this realparton was emitted.
   */
  set<pair<pair<tRealPartonPtr, tRealPartonPtr>,
	   pair<tRealPartonPtr, tRealPartonPtr> > > exchanges;
  set<pair<pair<tRealPartonPtr, tRealPartonPtr>,
	   pair<tRealPartonPtr, tRealPartonPtr> > > cExchanges;

  /**
   * if the parton is a valence parton in the last saved state, and the info.
   */
  bool cValence;
  Energy cValencePlus;
  TransverseMomentum cValencePT;

  /**
   * The parton that inherited the valence status from this parton when
   * it got set to NO. OR the parton that it inherited the valence status from.
   */
  tRealPartonPtr movedValence;
  tRealPartonPtr cMovedValence;

  /**
   * The p- that has been given to this parton from the other state.
   **/
  Energy givenMinus;
  Energy cGivenMinus;

  /**
   * If the parton has been checked for DGLAP. No point in checking twice.
   **/
  bool DGLAPchecked;

  /**
   * Adds the interacting dipole to this parton, and recursively to its oMothers.
   */
  void addInteraction(tDipolePtr);

  /**
   * Return the pT the gluon got when emitted from its parents.
   */
TransverseMomentum opT() ;

  /**
   * Controls if the parton is ordered with respect to its mothers.
   * Valence partons return true.
   */
  bool isOrdered();

  /**
   * Controls if the parton is ordered with respect to its first mothers.
   * Valence partons return true.
   */
  bool isFirstOrdered();

  /**
   * Controls if the parton is ordered with respect to its second mothers.
   * Valence partons return true.
   */
  bool isSecondOrdered();

  /**
   * Controls if the parton is ordered with respect to its mother.
   * Valence partons return true.
   */
  bool isSingleOrdered();

  /**
   * Controls that the emission had enough p+ from the mothers
   **/
  bool hasEnergy();

  /**
   * What pT scale the inconsistencies of the parton is associated with.
   * Returns ZERO if parton is ok.
   * problemScale return the max of the other scales.
   **/
  Energy problemScale();
  Energy orderedScale();
  Energy DGLAPScale();
  Energy motherScale();

  /**
   * Some of the scales above has to be calculated at emission, remember/set these here.
   **/
  Energy emissionScale;
  tRealPartonPtr emissionCause;
  void checkEmissionProblem();

  /**
   * Returns the cause of the problem with the parton.
   **/
  RealPartonPtr findCause();
  RealPartonPtr orderedCause();
  RealPartonPtr DGLAPCause();
  RealPartonPtr motherCause();

  /**
   *  Finds the main suspect for causing this realparton to be unordered.
   **/
  RealPartonPtr findSuspect();

  /**
   * Controls if the parton has any children.
   */
  bool hasChild() const {
    return !children.first.empty() || !children.second.empty();
  }

  /**
   * Controls if the parton should be absorb for non-DGLAP chain,
   * using both future and history of the parton.
   * Valence partons return true.
   */
  bool DGLAPSafe(InvEnergy scale = ZERO);
  bool singleDGLAPSafe(InvEnergy scale = ZERO);

  /**
   * Tries to anticipate if the parton will get DGLAPVetoed
   * using only the backwards history.
   * Does not give false positives.
   */
  bool willBeDGLAPVetoed();

  /**
   * Randomly
   */
  bool checkDGLAPSafe(tRealPartonPtr, tRealPartonPtr, InvEnergy scale = ZERO);

  /**
   * Removes the mother links and undoes the recoils.
   */
  void eraseMothers();

  /**
   * Removes the children links and undoes the recoils.
   */
  void eraseChildren();

  /**
   * Finds the colour neighbors.
   */
  RealPartonPtr firstColourNeighbor();
  RealPartonPtr secondColourNeighbor();

  /**
   * Finds and returns the (guessed) colour neighbours within range.
   * Looks among mothers and childs on first or second side.
   * Due to the swing, this may not reproduce the original neighbor.
   */
  RealPartonSet effectiveParton(InvEnergy range, bool firstSide);

  /**
   * The effective plus and minus of the particles within range.
   * @first flags to look for partons at first or second side.
   */
  pair<Energy, Energy> effectivePlusMinus(InvEnergy range, bool firstSide);

  /**
   * The amount of minus required to put all partons within range on shell.
   */
  Energy inRangeMinus(InvEnergy range, bool firstSide);

  /**
   * Searches through the effective parton and return true if it finds any negative p+ or p-.
   **/
  bool searchNegative(InvEnergy range, bool firstSide);

  /**
   * Does the recoil for the parton p. Uses different schemes from EventFiller::recoilScheme.
   */
  void doRecoil( RealPartonPtr p, Energy recPlus, TransverseMomentum recPT );

  /**
   * Recoils p with constant p+.
   */
  void doXmasRecoil( RealPartonPtr p, Energy recPlus, TransverseMomentum recPT );

  /**
   * Recoils p with constant y. Extra p+ taken from parents.
   */
  void doRainfallRecoil( RealPartonPtr p, Energy recPlus, TransverseMomentum recPT );

  /**
   * Recoils p straight towards (y=0, log(pT)=log(W/2)). Takes extra p+ from parents.
   **/
  void sunshineRecoil(RealPartonPtr p, Energy recPlus, TransverseMomentum recPT);

  /**
   * coherently recoils the parton p and the real partons within range. plus and recoil refers
   * to how much is taken away from the partons, not how much they recieve.
   * if not enough plus in the partons, y and p- are not updated.
   * firstSide indicates if the recoil is the first or second mother.
   */
  void doEffectiveRecoil( RealPartonPtr p, InvEnergy range, bool firstSide,
			  Energy recPlus, TransverseMomentum recPT );

  /**
   * Same as above, but forced to weight the partons after p+
   */
  void doPlusWeightedRecoil( RealPartonPtr p, InvEnergy range, bool firstSide,
			  Energy recPlus, TransverseMomentum recPT );

  /**
   * undoes all recoils made when this real parton was emitted.
   */
  void undoRecoils();

  /**
   * undoes all recoils this real parton caused in the colliding state.
   */
  void undoInteractionRecoils();

  /**
   * Sets NO and updates the relatives. returns false.
   */
  bool setNO();

  /**
   * Sets YES and updates the relatives. returns true if ordered.
   */
  bool setYES();

  /**
   * Sets YES and updates the relatives. returns true if enough energy.
   */
  bool quickSetYES();

  /**
   * Updates y and minus, decide from pT and plus.
   */
  void updateYMinus();

  /**
   * removes all valence substitutes' valence status, and sets this to valence.
   */
  void reclaimValenceStatus();

  /**
   * Undoes all recoilswings that were triggered when the real parton was emitted.
   */
  void eraseExchanges();

  /**
   * Finds the youngest child emitted before the argument was emitted.
   * If none, the corresponding mother is returned if checkMothers.
   */
  tRealPartonPtr youngestFirstChild(tRealPartonPtr, bool checkMother);
  tRealPartonPtr youngestSecondChild(tRealPartonPtr, bool checkMother);

  /**
   * Sets the mothers of the parton by searching for YES backwards.
   */
  void setMothers();

  /**
   * Sets the single mother of the parton by calling findmothers.
   */
  void setMother();

  /**
   * set the mothers, and change the kid of the mother.
   */
  void setFirstOMother(tRealPartonPtr);
  void setSecondOMother(tRealPartonPtr);
  void setOMother(tRealPartonPtr);

  /**
   * removes the given child, and undoes the recoil.
   */
  void eraseFirstChild(tRealPartonPtr);
  void eraseSecondChild(tRealPartonPtr);

  /**
   * Set the momentum using the valencePT and oY in theParton.
   * Should be the starting momentum.
   */
  void setValenceMomentum();

  /**
   * sets p_mu of this realparton, and recoils the mothers.
   * Return false if not 2 mothers that are ordered.
   */
  bool doRecoil();

  /**
   * sets p_mu of this realparton, and recoils the mother.
   * Returns false if the mother is not ordered.
   */
  bool doSingleRecoil();

  /**
   * sets p_mu of this realparton, and recoils the mothers.
   * Return false if not enough energy to do emission.
   */
  bool quickDoRecoil();

  /**
   * check if its mothers were swinged together, and in that maybe gives them an
   * extra recoil. returns true if recoil is done.
   */
  bool doSwingRecoil();

  /**
   * sets the parton, and recursively its parents, as onShell.
   */
  void setOnShell();

  /**
   * Moves the valence status of this parton to the first child.
   * returns false if no kids, or if kids are already valence partons.
   */
  bool moveValenceStatus();

  /**
   * Moves the valence status to the realparton rp.
   */
  bool moveValenceStatus(RealPartonPtr);

  /**
   * Returns the oldest child, or grandchild etc, that is younger than the argument.
   */
  tRealPartonPtr oldestHeir(tRealPartonPtr);

  /**
   * Returns a set of all ancestors, the mothers, the mothers mothers, etc.
   */
  RealPartonSet ancestors();

  /**
   * Checks how much minus is missing from the parton and its parents,
   * not checking the ones in @checked.
   */
  Energy checkMinus(RealPartonSet checked);

  /**
   * Checks how much minus is missing from the effective parton and its parents,
   * not checking the ones in @checked.
   */
  Energy effectiveCheckMinus(InvEnergy range, bool firstSide, RealPartonSet checked);

  /**
   * Gives the parton and it's mothers the minus it needs, and returns how much was missing.
   */
  Energy giveMinus();

  /**
   * Gives the effective parton and it's mothers the minus it needs, and returns how much was missing.
   */
  Energy effectiveGiveMinus(InvEnergy range, bool firstSide);

  /**
   * Emits a recoiler gluon with transverse momentum rec, taking a ratio plusRatio of the p+
   */
  void emitRecoiler(TransverseMomentum rec, double plusRatio);

  /**
   * Returns a list of relative weights for the partons in the effective parton. decides how much
   * recoil each parton should take.
   */
  list<double> effectiveWeights(const RealPartonSet & partons);

  /**
   * Saves current configuration.
   */
  void saveState();

  /**
   * reverts to last saved configuration.
   */
  void revert(tDipolePtr intDip);

  /**
   * Returns true if the pT fits with wath is expected from the relatives.
   */
  bool checkMomentum() const;

  /**
   * the rapidity at which the parton was originally emitted.
   */
  const double oY() const {
    return theParton->oY();
  }

  /**
   * For debugging, can all be removed
   */
  TransverseMomentum intRecoil;
  TransverseMomentum cIntRecoil;
  mutable int plus1veto;
  mutable int plus2veto;
  mutable int minus1veto;
  mutable int minus2veto;
  mutable int unDGLAP;
  void checkPS() const;



  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RealParton & operator=(const RealParton &);

};

}

#endif /* DIPSY_RealParton_H */
