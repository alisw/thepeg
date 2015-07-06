// -*- C++ -*-
#ifndef DIPSY_RealPartonState_H
#define DIPSY_RealPartonState_H
//
// This is the declaration of the RealPartonState class.
//

#include "ThePEG/Config/ThePEG.h"
#include "RealPartonState.fh"
#include "RealParton.h"
#include "Parton.h"
#include "Dipole.h"
#include "DipoleState.fh"

#include <iostream>

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the RealPartonState class.
 */
class RealPartonState: public Base {

public:

  /**
   * Get typedef from RealParton for convenience.
   */
  typedef RealParton::RealPartonSet RealPartonSet;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  RealPartonState();

  /**
   * The destructor.
   */
  virtual ~RealPartonState();
  //@}

public:

  /**
   * A map relating all generated partons to the corresponding
   * RealParton objects.
   */
  map<tPartonPtr,RealPartonPtr> partons;

  /**
   * The p+ and p- that the original particle had.
   */
  Energy plus;
  Energy minus;

  /**
   * The p- missing that is not reflected in the valence partons.
   */
  Energy minusDeficit;

  /**
   * These are the partons that need to be checked to find a
   * consistent evolution history,
   */
  RealPartonSet toCheck;

  /**
   * The prime suspect of making the state not consistent.
   **/
  RealPartonSet suspects;

  /**
   * These partons are considered to be valence and do not need to be
   * checked.
   */
  RealPartonSet valence;
  RealPartonSet oValence;

  /**
   * These are the partons which would be directly involved in the
   * interaction.
   */
  RealPartonSet interacting;

  /**
   * These are the dipoles which would interact if we can find a
   * consistent evolution history.
   */
  list<tDipolePtr> interactions;
  list<pair<bool, bool> > doesInts;

  /**
   * Tested interacting dipoles that failed as primary interactions.
   */
  set<tDipolePtr> failedInts;

  /**
   * True if a consistent evolution history has been found for this
   * state.
   */
  bool consistent;

  /**
   * True if the last treated interaction was with a previously interacting dipole.
   */
  bool rescatter;

  /**
   * The set of real partons that recoiled over y0 in the last interaction.
   */
  RealPartonSet backwardsPartons;

  /**
   * The set of virtual fluctuations in the state.
   **/
  vector<RealPartonSet> flucts;
  vector<RealPartonSet> cFlucts;

  /**
   * Adds an interaction to the state, updating all partons and toCheck.
   */
  void newInteraction(tDipolePtr intDip, tDipolePtr otherIntDip, bool firstInt, bool secondInt,
		      Energy rec1 = ZERO, Energy rec2 = ZERO);

  /**
   * Adds a real parton to the state, and to toCheck.
   * Recurs if it has no previous interactions.
   */
  tRealPartonPtr addParton(tPartonPtr);

  /**
   * Adds all valence partons in the state to the real state.
   */
  void addValence(DipoleState & state);

  /**
   * Finds the RealParton associated with p, or creates one if it doesnt exist.
   */
  tRealPartonPtr getReal(tPartonPtr p);

  /**
   * Assuming the past fix, looks for a consistent forward history in which the
   * parton can exist.
   */
  bool findConsistentEvolution(RealPartonSet::iterator);

  /**
   * Does the evolution of toCheck wihtout checking ordering etc.
   * Only stops if a parton doesnt have enough p+, otherwise goes all the way through.
   **/
  void doEvolution();

  /**
   * checks a parton for non-DGLAP and fixes by merging if needed.
   * this method recurs backwards to check the pt max first.
   * Compares to the scale provided, or find it itself if zero.
   **/
  bool checkFixDGLAP(RealPartonPtr rp, InvEnergy scale = ZERO);
  /**
   * as above, but for partons with a single mother
   **/
  bool checkFixSingleDGLAP(RealPartonPtr rp, InvEnergy scale = ZERO);

  /**
   * checks if mom has a DGLAP-checkable structure, and the mother scale
   * corresponds to a pt larger than the provided scale. returns false
   * if either of those fails.
   **/
  bool isMotherStepUp(tRealPartonPtr mom, bool firstSide, InvEnergy scale);

  /**
   * merges an FSR double counted parton with its mother.
   **/
  void fixFSRDoubleCounting(tRealPartonPtr rp);

  /**
   * checks if the parton couldve been emitted as FSR from its mother and children
   * does not check if the child is interacting.
   **/
  bool inFSRRegion(tRealPartonPtr rp);

  /**
   * takes plus from the parents of rp to set give rp a plus of zero.
   * does nothing if rp already has positive plus
   **/
  Energy fixNegativePlus(RealPartonPtr rp);

  /**
   * checks if rp has a mother in one colour direction, and a child in the other.
   * this means that the parton should be checked for DGLAP.
   **/
  bool isOutside(tRealPartonPtr rp, bool firstSide);

  /**
   * The distance to the last child (or int). Recurs if merged.
   **/
  InvEnergy childScale(tRealPartonPtr rp, bool firstSide);

  /**
   * Fixes a DGLAP inconsistency by merging.
   **/
  void fixDGLAP(RealPartonPtr rp);

  /**
   * Fixes a DGLAP inconsistency by merging in the singlemother mode.
   **/
  void fixSingleDGLAP(RealPartonPtr rp);

  /**
   * Fixes a DGLAP inconsistency by merging.
   **/
  void fixUnOrdered(RealPartonPtr rp, bool forced = false);

  /**
   * Adds all plus and pT in each virtual set to the first parton,
   * and flags the rest as off shell. Manually sets NO without calling setNO(),
   * to avoid other recoils being canceled outside the virtual set.
   * should be called only last, after all evos and interactions are done.
   **/
  void mergeVirtuals();

  /**
   * changes all recoils to the "cloud scheme", where the recoil kicks out a fraction
   * of the parent in a new recoiler gluon.
   **/
  void addRecoilers();

  /**
   * changes rp to the "cloud scheme", where recoils kick out a fraction
   * of the parent in a new recoiler gluon.
   **/
  void addRecoiler(RealPartonPtr rp);
  /**
   * a different (simpler) implementation of the same thing.
   **/
  void addRecoiler2(RealPartonPtr rp);

  /**
   * checks if there are any YES partons with negative plus or minus.
   **/
  bool checkForNegatives();

  /**
   * checks if all fluctuations contains partons with identical momenta.
   * (for debugging)
   **/
  bool checkFlucts();

  /**
   * couts the total and onshell plus and minus.
   * (for debugging)
   **/
  bool checkPlusMinus();

  /**
   * Joins the partons in a virtual pair, sharing their momenta and flaging them
   * for reabsorbtion after interaction.
   **/
  void merge(RealPartonPtr rp1, RealPartonPtr rp2);

  /**
   * Splits the fluctuation so that rp1 and rp2 will not be merged.
   * Does nothing if they are not in a common fluctuation.
   **/
  void splitFluct(RealPartonPtr rp1, RealPartonPtr rp2);

  /**
   * Divides up the plus and pT evenly and then updates minus.
   **/
  void makeCollinear(const RealPartonSet & rps);

  /**
   * Divides up the plus weighted by exp(-oY),
   * with oY the original y the parton was emitted at.
   **/
  void redistributePlus(const RealPartonSet & rps);

  /**
   * Calls setNO for all partons, except the original valence partons, that get setYES.
   */
  void reset();

  /**
   * Calls setNO for all non-valence in toCheck. Calls setYES for valence in toCheck.
   */
  void undoEmissions();

  /**
   * Sets all YES partons in toCheck after rp in oY to NO, and undoes recoils by YES or INT.
   */
  void setFutureNO(tRealPartonPtr rp);

  /**
   * Sets all YES partons in rps to NO, calling the setNO of each parton in reverse oY order.
   */
  void setNO(const RealPartonSet & rps);

  /**
   * Looks through the partons in toCheck for the problem with highest pT scale.
   **/
  RealPartonPtr findWorstProblem();

  /**
   * Tries to solve the provided problem by removing a parton.
   * Returns false if no solution is found.
   **/
  bool fix(RealPartonPtr problem);

  /**
   * Adds the interacting dipole and checks if the evolution is self consistent.
   */
  bool controlEvolution(tDipolePtr intDip, tDipolePtr otherIntDip);

  /**
   * Adds the interacting dipole and does a quick check to remove unordered partons.
   */
  bool singleControlEvolution(tDipolePtr intDip, tDipolePtr otherIntDip,
			      bool firstInt, bool secondInt, Energy rec1 = ZERO, Energy rec2 = ZERO);

  /**
   * Finds an evolution by eliminating one parton at a time.
   */
  bool nonRecursiveControlEvolution(tDipolePtr intDip, tDipolePtr otherIntDip);

  /**
   * Adds the interacting dipole and remakes the consistency check
   * from the start.
   */
  bool fullControlEvolution(tDipolePtr intDip, tDipolePtr otherIntDip,
			    bool firstInt, bool secondInt, Energy rec1 = ZERO, Energy rec2 = ZERO);

  /**
   * adds all partons which have at least one inteaction to the partons to be checked.
   */
  void checkAllInteracting();

  /**
   * adds partons that belong to only the provided interaction to toCheck.
   */
  void checkOnlyNew(tDipolePtr intDip);

  /**
   * adds partons that belong to the provided interaction to toCheck.
   */
  void checkInteraction(tDipolePtr intDip);

  /**
   * adds rp and its oMothers (recursively) to toCheck.
   */
  void checkHistory(tRealPartonPtr rp);

  /**
   * adds the future of every realparton in toCheck to toCheck.
   */
  void checkFuture();

  /**
   * Tells all partons in toCheck to save away current status as
   * consistent.
   */
  void saveState();

  /**
   * Reverts the state and all its parton to last saved state, and removes last interaction
   */
  void revertToPrevious(DipolePtr intDip);

  /**
   * Removes the correct interaction tags on the partons after a failed rescattering.
   */
  void removeLastRescatter(tDipolePtr intDip);

  /**
   * Removes/adds the interaction tags on the parton and all its ancestors.
   */
  void removeInt(tRealPartonPtr rp, DipolePtr intDip);
  void addInt(tRealPartonPtr rp, DipolePtr intDip);


  // /**
  //  * Return a list of how much p- each interaction needs to put the state on shell.
  //  * NOT USED??
  //  */
  // list<Energy> neededMinuses();

  /**
   * Checks how much p- all valence partons need.
   **/
  Energy neededValenceMinus();

  /**
   * Sets onShell to true for all YES partons deriving from this interaction.
   */
  void setOnShell(tDipolePtr intDip);

  /**
   * Cleans up as if the interaction did not happen by setting both
   * partons NO.  Assumes no rescattering.
   */
  void removeInteraction(tDipolePtr);

  /**
   * Removes the partons in backwardsPartons and finds a new parton to do the interaction.
   */
  void removeBackwards();

  /**
   * Prints the state to "realState.dat", and can pauses the program.
   */
  void plotState(bool pause) const;

  /**
   * Prints this state and rrs (mirrored) to "realState.dat", and can pauses the program.
   */
  void plotBothStates(RealPartonStatePtr rrs, bool pause) const;

  /**
   * Prints the state to @file, and can pauses the program.
   * The file has to be opened and closed outside this function.
   * The state can be mirrored in y = 0 in the plot.
   */
  void plotStateToFile(ostream & file, bool mirror) const;

  /**
   * returns thenumber of YES and NO particles in the state.
   */
  pair<int, int> countYESNO() const;

  /**
   * Returns the average length in rapidity between a parton and its mother.
   */
  double avYInEvo() const;

  /**
   * Runs some tests on the state to check for inconsistencies.
   */
  bool diagnosis(bool pause) const;

  /**
   * for debugging.
   */
  mutable int calls;
  mutable int unOrdered;
  mutable int unDGLAP;
  mutable double currentY;
  mutable int partonCalls;
  mutable bool monitored;
  mutable TransverseMomentum totalRecoil;
  mutable TransverseMomentum cTotalRecoil;

  Energy totalMinus();
  Energy totalPlus();

protected:

  /**
   * Exception class to signal bad kinematics.
   */
  struct RealPartonKinematicException: public Exception {};

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RealPartonState & operator=(const RealPartonState &);

};

}


#endif /* DIPSY_RealPartonState_H */
