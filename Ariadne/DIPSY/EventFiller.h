// -*- C++ -*-
#ifndef DIPSY_EventFiller_H
#define DIPSY_EventFiller_H
//
// This is the declaration of the EventFiller class.
//

#include "EventFiller.fh"
#include "ThePEG/Handlers/HandlerBase.h"
#include "DipoleEventHandler.fh"
#include "RealPartonState.fh"
#include "RealParton.fh"
#include "DipoleState.h"
#include "ImpactParameters.h"
#include "DipoleXSec.h"
#include "DipoleAbsorber.h"
#include "ThePEG/Utilities/Current.h"


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
 * The EventFiller class is able to produce an initial
 * ThePEG::Collision from two colliding DipoleStates.
 *
 * @see \ref EventFillerInterfaces "The interfaces"
 * defined for EventFiller.
 */
class EventFiller: public HandlerBase {

public:

  /**
   * Copy FList typedef from DipoleXSec.
   */
  typedef DipoleXSec::FList FList;

  /**
   * A String is simply a vector of colour-connected partons.
   */
  typedef DipoleState::String String;

  /**
   * A vector of dipole pairs represented as iterators into an FList.
   */
  typedef DipoleXSec::DipolePairVector DipolePairVector;

  /**
   * An ordered map of dipole pairs represented as iterators into an
   * FList.
   */
  typedef DipoleXSec::DipolePairMap DipolePairMap;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  EventFiller();

  /**
   * The destructor.
   */
  virtual ~EventFiller();
  //@}

public:

  /** @name Virtual functions which can be overridden in subclasses. */
  //@{
  /**
   * Fill the current collision object in the given DipoleEventHandler
   * with the final state gluons produced when dipole states \a dl and
   * \a dr collides with impact parameters \a b.
   * @return the total interaction probability.
   */
  virtual double fill(Step & step, DipoleEventHandler & eh, tPPair inc,
		      DipoleState & dl, DipoleState & dr,
		      const ImpactParameters & b) const;

  /**
   * Select which dipole pairs should interact, and return them in the
   * order in which they should be processed.
   */
  virtual pair<RealPartonStatePtr, RealPartonStatePtr>
  selectInteractions(const FList & fl, const ImpactParameters & b, const DipoleXSec & xSec) const;



 
  /**
   * Tries to do the interaction between the real states lrs and rrs with the two dipoles
   * pointed to by inter, assuming that the interactions in inters already has been tested
   * and accepted. States collide at impact paramter b and using the recoils in xSec.
   */
  bool addInteraction(FList::const_iterator inter,
		      RealPartonStatePtr lrs, RealPartonStatePtr rrs,
		      DipolePairVector & inters,
		      const ImpactParameters & b, const DipoleXSec & xSec) const ;

  /**
   * Extract all strings after according to the given interactions.
   */
  virtual vector<String>
  extractStrings(DipoleState & dl, DipoleState & dr,
		 pair<RealPartonStatePtr, RealPartonStatePtr>,
		 const ImpactParameters & b) const;

  /**
   * Fill the given Step with a SubProcess object using the given
   * incoming particles and list of strings.
   */
  virtual bool fillStep(Step & step, tPPair incoming,
			const vector<String> & strings) const;

  /**
   * Fix up valens partons if they were not of the correct flavour or
   * if they should collapse into a hadron.
   */
  //  virtual void fixValence(Step & step, DipoleState & dl, DipoleState & dr) const;

  /**
   * get the recoil scheme.
   */
  inline int recoilScheme() const {
    return theRecoilScheme;
  }

  /**
   * set the recoil scheme.
   */
  inline void recoilScheme(int x) {
    theRecoilScheme = x;
  }

  /**
   * get the mode.
   */
  inline int mode() const {
    return theMode;
  }

  /**
   * set the mode.
   */
  inline void mode(int x) {
    theMode = x;
  }

  /**
   * get the singleMother.
   */
  inline int singleMother() const {
    return theSingleMother;
  }

  /**
   * set the singleMother.
   */
  inline void singleMother(int x) {
    theSingleMother = x;
  }

  /**
   * get the DGLAPinPT.
   */
  inline int DGLAPinPT() const {
    return theDGLAPinPT;
  }

  /**
   * set the DGLAPinPT.
   */
  inline void DGLAPinPT(int x) {
    theDGLAPinPT = x;
  }
  
  /**
   * get the ValenceChargeNormalisation
   **/
  inline int valenceChargeNormalisation() const {
    return theValenceChargeNormalisation;
  }

  /**
   * The minimum squared invariant transverse momentum allowed for a
   * gluon in a string.
   */
  Energy2 pT2Cut() const {
    return sqr(thePTCut);
  }

  /**
   * Return the squared invariant transverse momentum of gluon \a i2
   * in a string given its colour-neighbours \a i1 and \a i3.
   */
  static Energy2 invPT2(const String & str, int i1, int i2, int i3);

  /**
   * Remove a gluon \a i2 in a string shuffling its momenum to its
   * colour-neighbours \a i1 and \a i3.
   */
  static bool removeGluon(const String & str, int i1, int i2, int i3);

  /**
   * get the effectiveWeights.
   */
  inline int effectiveWeights() const {
    return theEffectiveWeights;
  }


  /**
   * get the FS swing time.
   */
  inline double FSSwingTime() const {
    return theFSSwingTime;
  }

  /**
   * get the step size for the FS swing time.
   */
  inline double FSSwingTimeStep() const {
    return theFSSwingTimeStep;
  }

  //@}

public:

  /**
   * Get the object used to absorb non-interacting dipoles.
   */
  inline DipoleAbsorberPtr absorber() const {
    return theAbsorber;
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

  /**
   * The weight which will be used for the event being generated.
   */
  mutable double currentWeight;
  mutable bool fail;

  /**
   * Debug/testing
   */
protected:
  mutable double nInt;
  mutable double nInt1;
  mutable double nTried;
  mutable double nEvents;
  mutable double nTestedInts;
  mutable double foundInt;
  mutable double failedEvo;
  mutable double failedRec;
  mutable double rescatter;
  mutable double overTenEvents;
  mutable double rejectedEvents;
  mutable double nLoops;
  mutable double nYes;
  mutable double nNo;
  mutable double avYInEvo;
  mutable double avYInInt;

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


protected:

  /**
   * balances p+/p- by moving the entire state in rapidity. Assumed to be
   * UNMIRRORED states, ie, plus and minus are flipped for the right state.
   **/
  void fixBoost(RealPartonStatePtr lrs, RealPartonStatePtr rrs) const;

  /**
   * Checks through the final state and tries to fix the most blatant errors.
   * Specially tries to catch bad momenta that can cause crashes.
   **/
  void dodgeErrors(DipoleStatePtr finalState) const;

  /**
   * Checks if the interacting partons in the current states have enough 
   * p+- to set the other state on shell.
   */
  bool controlRecoils(DipolePairVector & sel, 
		    RealPartonStatePtr rrs, RealPartonStatePtr lrs,
		      const ImpactParameters & b, const DipoleXSec & xSec,
		      pair<pair<bool, bool>, pair<bool, bool> > doesInt) const;

  /**
   * Removes the off shell partons and recouples colour flow. p_mu unchanged.
   */
  void removeVirtuals(DipoleStatePtr state) const;

  /**
   * Removes the parton and recouples colour flow. p_mu unchanegd.
   */
  void removeParton(tPartonPtr p) const;

  /**
   * find the pT scale through the dipolestates eventhandlers emitter object.
   */
  double pTScale(DipoleState &) const;

  /**
   * Takes statistics on the number of participants
   * in a heavy ion collision.
   */
  void countParticipants(const DipoleState & dl, const DipoleState & dr, const InvEnergy b) const;

  /**
   * The object used to absorb non-interacting dipoles.
   */
  DipoleAbsorberPtr theAbsorber;

  /**
   * What scheme to use when doing recoils.
   */
  int theRecoilScheme;

  /**
   * In what mode to create the real state. Fast vs consistent.
   */
  int theMode;

  /**
   * If partons have one or two mothers.
   */
  int theSingleMother;

  /**
   * If DGLAP supression is made in terms of pt rather than r.
   */
  int theDGLAPinPT;

  /**
   * How to distribute recoil among the members of an effective parton.
   */
  int theEffectiveWeights;

  /**
   * How long time (in GeV-1) the FS get's to do colour reconnections.
   */
  double theFSSwingTime;

  /**
   * How large time steps (in GeV-1) used in the FS colour
   * reconnections.
   */
  double theFSSwingTimeStep;

  /**
   * How the valence charge is handled.
   **/
  int theValenceChargeNormalisation;

  /**
   * The minimum invariant transverse momentum allowed for a gluon in a string.
   */
  Energy thePTCut;

  /**
   * Options for removing gluons with too small invariant transverse
   * momentum.
   */
  int theSoftRemove;

public:

  /**
   * Exception class for space-like gluon momenta.
   */
  struct SpaceLikeGluons: public Exception {};

  /**
   * Exception class for failed gluon removal.
   */
  struct RemoveGluonException: public Exception {};

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EventFiller & operator=(const EventFiller &);

};

}

#endif /* DIPSY_EventFiller_H */
