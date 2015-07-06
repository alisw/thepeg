// -*- C++ -*-
#ifndef DIPSY_DipoleXSec_H
#define DIPSY_DipoleXSec_H
//
// This is the declaration of the DipoleXSec class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "DipoleXSec.fh"
#include "Dipole.fh"
#include "DipoleState.fh"
#include "ImpactParameters.h"
#include "RealPartonState.fh"
#include "RealParton.fh"
#include "EffectiveParton.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * DipoleXSec is the base class of all objects capable of calculating
 * the scattering probability for two colliding dipoles. The main
 * virtual functions to be overridden are called fij, sumf and flist.
 *
 * @see \ref DipoleXSecInterfaces "The interfaces"
 * defined for DipoleXSec.
 */
class DipoleXSec: public HandlerBase {

public:

  /**
   * A maptype used to store an ordered list of scattering
   * probabilities (the bare one paired with the unitarized one)
   * between dipoles in colliding dipole systems.
   */
  typedef multimap<pair<double,double>, pair<tDipolePtr,tDipolePtr>,
                   std::greater< pair<double,double> > > FList;

  /**
   * Interaction pT of the four partons involved in an interaction. Default order is
   * the first and second parton in colour order from the left state first,
   * then the first and second parton from right state.
   */
  typedef pair< pair<TransverseMomentum, TransverseMomentum>,
		pair<TransverseMomentum, TransverseMomentum> > InteractionRecoil;

  /**
   * A vector of dipole pairs represented as iterators into an FList.
   */
  typedef vector<FList::const_iterator> DipolePairVector;

  /**
   * An ordered map of dipole pairs represented as iterators into an
   * FList.
   */
  typedef multimap<double,FList::const_iterator,std::greater<double> >
          DipolePairMap;

  /**
   * A interaction of RealPartons, with all needed information about the interaction.
   */
  struct RealInteraction {
    RealPartonStatePtr lrs, rrs;
    DipolePtr d1, d2;
    RealPartonPtr p11, p12, p21, p22;
    InvEnergy range11, range12, range21, range22;
    Energy effectivePlus11, effectivePlus12, effectivePlus21, effectivePlus22;
    Energy effectiveMinus11, effectiveMinus12, effectiveMinus21, effectiveMinus22;
    bool max11, max12, max21, max22;
    double P11, P12, P21, P22;
  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline DipoleXSec();

  /**
   * The copy constructor.
   */
  inline DipoleXSec(const DipoleXSec &);

  /**
   * The destructor.
   */
  virtual ~DipoleXSec();
  //@}

public:

  /** @name Virtual functions to be overridden by subclasses. */
  //@{
  /**
   * Calculate the scattering probability for the two given dipoles
   * using the given ImpactParameters.
   */
  virtual double fij(const pair<tPartonPtr, tPartonPtr>,
		     const pair<tPartonPtr, tPartonPtr>,
		     const ImpactParameters &  b,
		     bool veto = true) const;

  /**
   * Calculate the scattering probability for the two given dipoles
   * using the given ImpactParameters.
   */
  virtual double fij(const Dipole &,
		     const Dipole &,
		     const ImpactParameters &  b,
		     bool veto = true) const;

  /**
   * Return false if an interaction would be kinematically forbidden.
   */
  virtual bool kinematicsVeto(const pair<tPartonPtr, tPartonPtr>,
			      const pair<tPartonPtr, tPartonPtr>,
			      const ImpactParameters &  b) const;

  /**
   * Return false if an interaction would be kinematically forbidden.
   */
  virtual bool kinematicsVeto(const Dipole &,
			      const Dipole &,
			      const ImpactParameters &  b,
			      const pair<bool,bool> & ints) const;

  /**
   * Calculate the total scattering probability for the two given
   * dipole systems using the given ImpactParameters.
   */
  virtual double sumf(const DipoleState &, const DipoleState &,
		      const ImpactParameters &) const;

  /**
   * Calculate the total scattering probability all dipole pairs 
   * the two given dipole systems using the given ImpactParameters.
   */
  virtual FList flist(const DipoleState &, const DipoleState &,
		      const ImpactParameters &) const;

  /**
   * Calculate the total scattering probability all dipole pairs 
   * the two given dipole systems using the given ImpactParameters.
   * This version backtrack vetoed dipoles pairs and checks the parents.
   */
  virtual FList effectiveFlist(const DipoleState &, const DipoleState &,
			       const ImpactParameters &) const;

  /**
   * Return a unitarized scattering probability, given the
   * ununitarized one.
   */
  virtual double unitarize(double f) const;

  /**
   * Returns the recoils the interaction would like to give the 4 involved partons.
   * Assumes the states are not yet rotated.
   */
  virtual InteractionRecoil
  recoil(const pair<tPartonPtr, tPartonPtr> left,
	 const pair<tPartonPtr, tPartonPtr> right,
	 const ImpactParameters & b,
	 pair<pair<bool, bool>, pair<bool, bool> > doesInt,
	 pair<bool,bool> ints0) const;

  /**
   * Returns the recoils the interaction would like to give the 4 involved partons.
   * Assumes the states are not yet rotated.
   */
  virtual InteractionRecoil
  recoil(const pair<tPartonPtr, tPartonPtr> left,
	 const pair<tPartonPtr, tPartonPtr> right,
	 const ImpactParameters & b,
	 pair<pair<bool, bool>, pair<bool, bool> > doesInt = 
	 make_pair(make_pair(true, true), make_pair(true, true))) const;

  /**
   * Does the transverse recoil on the four partons.
   */
void doTransverseRecoils(RealInteraction i, InteractionRecoil recs) const;

  /**
   * creates and initialises a RealInteraction with realpartons, ranges and max.
   */
  virtual RealInteraction
  initialiseInteraction(const pair<DipolePtr, DipolePtr> inter,
			RealPartonStatePtr lrs, RealPartonStatePtr rrs,
			pair<pair<bool, bool>, pair<bool, bool> > doesInt,
			const ImpactParameters & b) const;

  /**
   * Decides which of the four partons actually interact.
   */
  virtual pair<pair<bool, bool>, pair<bool, bool> >
  doesInt(const pair<tPartonPtr, tPartonPtr> left,
	  const pair<tPartonPtr, tPartonPtr> right,
	  const ImpactParameters & b) const;

  /**
   * Selects the two partons (out of the four) that are used in theInteraction == 0.
   */
  virtual pair<bool, bool> int0Partons(tcPartonPtr p11, tcPartonPtr p12,
				       tcPartonPtr p21, tcPartonPtr p22, 
				       const ImpactParameters & b) const;

  /**
   * calculates and sets the effective plus and minus of the partons of the interaction
   **/
  virtual void updateMomenta(RealInteraction * interaction) const;

  // /**
  //  * undoes the last recoil for each of the four partons.
  //  */
  // virtual void undoLastRecoil(RealPartonPtr p11, RealPartonPtr p12,
  // 			      RealPartonPtr p21, RealPartonPtr p22) const;

  /**
   * Does the transverse recoil on the real partons with the pT supplied in recs.
   * Does the p+ alt. p- transfer of neededPlus nad neededMinus to put the
   * other state on shell. If not consistent, return false.
   * Also updates p_mu for the 4 effective partons. Assumes all are rightmoving.
   */
  virtual bool doInteraction(InteractionRecoil recs, const FList::const_iterator inter,
			     RealPartonStatePtr lrs, RealPartonStatePtr rrs,
			     pair<pair<bool, bool>, pair<bool, bool> > doesInt,
			     const ImpactParameters & b) const;

  /**
   * Tests the provided interactions in the same way doInteraction does, but
   * without transfering any momenta between the states.
   */
  virtual bool checkInteractions(DipolePairVector & ints, RealPartonStatePtr lrs,
				 RealPartonStatePtr rrs, const ImpactParameters & b) const;

  /**
   * undoes the recoils to the other state stored in all partons involved in the
   * interaction. Note that the p+- boosts are stored as recoils against themselves, and 
   * does not trigger here.
   **/
  virtual void undoPTRecoils(RealInteraction RI) const;

  /**
   * Performs the interactions provided. Returns false if something goes wrong.
   */
  virtual bool performInteractions(DipolePairVector & ints, RealPartonStatePtr lrs,
				   RealPartonStatePtr rrs, const ImpactParameters & b) const;

  /**
   * checks if the four interacting partns are ordered.
   * It is assumed that p2 not yet is mirrored in y=0.
   * The booleans indicate if the parton pairs are local max in pT.
   **/
  virtual bool checkOrderedInteraction(tRealPartonPtr p11, tRealPartonPtr p12,
				       tRealPartonPtr p21, tRealPartonPtr p22) const;

  /**
   * reconnects the colour flow in an interaction, taking care that rescatterings are
   * treated properly.
   */
  virtual bool reconnect(tDipolePtr d1, tDipolePtr d2) const;

  /**
   * checks the interactions in the two real states, ans picks out the ones that
   * did non-colour singlet exchanges.
   */
  virtual vector<pair<DipolePtr, DipolePtr> >
  getColourExchanges(tRealPartonStatePtr lrs, tRealPartonStatePtr rrs) const;

  /**
   * Return the interaction
   **/
  inline int interaction() const {
    return theInteraction;
  }

  //@}

  /**
   * The confinement scale.
   */
  InvEnergy rMax() const;

  /**
   * for debugging
   */
  mutable int Nfij;
  mutable int NScalVeto;
  mutable int NBVeto;
  mutable double scalVeto;
  mutable double bVeto;

private:

  /**
   * return the fractions of p+ and p- to be transfered to put everything on shell.
   */
  pair<double, double> findBoosts(Energy intPlus1, Energy intPlus2,
				  Energy intMinus1, Energy intMinus2,
				  Energy evoPlus2, Energy evoMinus1) const;

  /**
   * return the fractions of p+ and p- to be transfered to put everything on shell.
   **/
  pair<double, double> findBoosts(RealInteraction i) const;

  /**
   * boosts the interacting particles with the scales provided
   **/
  void doBoosts(RealInteraction i, pair<double, double> boosts) const;

  /**
   * boosts the affective partons with the scale x.
   */
  void doBoost(tRealPartonPtr p1, InvEnergy range1,
	       tRealPartonPtr p2, InvEnergy range2, double x) const;

  /**
   * reduce the recoil until the interacting partons are ordered in rapidity.
   **/
  void reduceRecoil(RealInteraction RI, InteractionRecoil recs) const;

  /**
   * checks if the interaction is ordered or not
   **/
  bool ordered(RealInteraction i, InteractionRecoil recs, const ImpactParameters & b) const;

  /**
   * sets all particles in the backwards cone on shell
   **/
  void setOnShell(RealInteraction i) const;

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

protected:

  /**
   * Exception class to signal bad kinematics.
   */
  struct InteractionKinematicException: public Exception {};

private:

  /**
   * The confinement scale.
   */
  InvEnergy theRMax;

  /**
   * Flag determining which interaction to be used.
   */
  int theInteraction;

  /**
   * What kind of kinematical ordering is requeired in the interaction.
   */
  int theIntOrdering;

  /**
   * What to do with large recoils, if anything.
   */
  int theRecoilReduction;

  /**
   * Check for off-shell incoming particles.
   */
  bool checkOffShell;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleXSec & operator=(const DipoleXSec &);

};

}

#include "DipoleXSec.icc"

#endif /* DIPSY_DipoleXSec_H */
