// -*- C++ -*-
#ifndef DIPSY_DipoleXSec_H
#define DIPSY_DipoleXSec_H
//
// This is the declaration of the DipoleXSec class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "DipoleXSec.fh"
#include "Dipole.h"
#include "DipoleState.fh"
#include "ImpactParameters.h"
#include "RealPartonState.fh"
#include "RealParton.fh"
#include "EffectiveParton.h"
#include "ShadowParton.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * A dipole-dipole interaction used for book keeping.
 */
struct DipoleInteraction {

  /**
   * Enumerate the resons for an interaction not working.
   */
  enum Status {
    UNKNOWN = -1, /**< Not yet checked. */
    ACCEPTED = 0, /**< Interaction is OK */
    PROPFAIL = 1, /**< Kinematics of incoming propagators failed. */
    KINEFAIL = 2, /**< Kinematics of interaction failed. */
    ORDERING = 3 /**< Ordering of interation failed. */
  };

  /**
   * Constructor taking two dipoles.
   */
  DipoleInteraction(Dipole & dlin, Dipole & drin,
		    const ImpactParameters &  bin, int ordering);

  /**
   * Prepare an interactions to be checked for the first time.
   */
  void prepare() const;

  /**
   * Check that the interaction can be performed. If \a mode >= 0 also
   * flag shadow partons on-shell. If \a mode > 0 also propagate
   * momenta to original partons.
   */
  Status check(int mode) const;

  /**
   * Accept the interaction.
   */
  void accept() const;

  /**
   * Reject this interaction.
   */
  void reject() const;

  /**
   * The two dipoles.
   */
  pair<tDipolePtr,tDipolePtr> dips;

  /**
   * The two dipoles next to these.
   */
  pair<tDipolePtr,tDipolePtr> dnext;

  /**
   * The two dipoles previous to these.
   */
  pair<tDipolePtr,tDipolePtr> dprev;

  /**
   * A pointer to the impact parameter object.
   */
  const ImpactParameters * b;

  /**
   * Which partons are actually interacting?
   */
  pair<tPartonPtr,tPartonPtr> ints;

  /**
   * Which partons are merely spectators
   */
  pair<tPartonPtr,tPartonPtr> spec;

  /**
   * Which partons are actually interacting?
   */
  mutable pair<tSPartonPtr,tSPartonPtr> sints;

  /**
   * The distance between the interacting partons.
   */
  InvEnergy2 d2;

  /**
   * The interaction strength (times two for total xsec).
   */
  double f2;

  /**
   * The unitarized interaction strength.
   */
  double uf2;

  /**
   * The transverse momentum recoil associated with the interaction.
   */
  TransverseMomentum rec;

  /**
   * Set to true if the two partons which interacts has already received a recoil.
   */
  mutable bool norec;

  /**
   * The transverse momentum scale associated with the interactions.
   */
  Energy kt;

  /**
   * The number of this interaction in the order of tried interactions.
   */
  mutable int id;

  /**
   * Current status of this interaction.
   */
  mutable Status status;

  /**
   * Ordering scheme for interactions
   */
  int intOrdering;

  /**
   * Check if the given parton is ordered wrt. the given light-cone
   * momenta.
   */
  bool disorder(bool doit, tcSPartonPtr p, Energy plus, Energy minus) const {
    return ( doit && ( p->plus() < plus || p->minus() > minus ) );
  }
  bool disorder(bool doit, Energy plus, Energy minus, tcSPartonPtr p) const {
    return ( doit && ( plus < p->plus() || minus > p->minus() ) );
  }

  /**
   * Check that the two scattered partons are ordered.
   */
  bool disorder(bool doit) const {
    return disorder(doit, sints.first,
		    sints.second->minus(), sints.second->plus());
  }

  /**
   * Return true if the given scales are ascending of decanding.
   */
  bool cending(Energy2 ptl, Energy2 ptm, Energy2 ptr) const {
    return ( ( ptl > ptm && ptm > ptr ) || ( ptl < ptm && ptm < ptr ));
  }

  /**
   * Return true if the interaction fails any of the ordering requirements.
   */
  bool orderfail(const ShadowParton::Propagator & ppl,
		 const ShadowParton::Propagator & ppr) const;

  void checkShadowMomentum(const LorentzMomentum & pin = LorentzMomentum()) const;

  /**
   * Print out sum of on-shell momenta.
   */
  void debug() const;

  /** Struct for ordering by dipoles. */
  struct DipoleOrder {
    /** Main operator */
    bool operator()(const DipoleInteraction & i1,
		    const DipoleInteraction & i2) const {
      return i1.dips.first < i2.dips.first ||
	( i1.dips.first == i2.dips.first &&
	  i1.dips.second < i2.dips.second );
    }
  };

  /** Struct for ordering by dipoles. */
  struct StrengthOrder {
    /** Main operator */
    bool operator()(const DipoleInteraction & i1,
		    const DipoleInteraction & i2) const {
      return i1.f2 > i2.f2;
    }
  };

  /**
   * A settype used to store an ordered list of potential interactions
   * between dipoles in colliding dipole systems. Ordered in
   * decreasing interaction strength.
   */
  typedef multiset<DipoleInteraction,DipoleInteraction::StrengthOrder> List;
  
  /** Struct for ordering by scale. */
  struct PTOrder {
    /** Main operator */
    bool operator()(const List::const_iterator & i1,
		    const List::const_iterator & i2) const {
      return i1->kt > i2->kt;
    }
  };

  /**
   * A settype used to store an ordered list of potential interactions
   * between dipoles in colliding dipole systems. Ordered in
   * decreasing interaction strength.
   */
  typedef multiset<List::const_iterator,PTOrder> PTSet;

  /**
   * A settype used to store a list of potential interactions
   * between dipoles in colliding dipole systems.
   */
  typedef set<DipoleInteraction,DipoleInteraction::DipoleOrder> DipList;

  /**
   * Check which orderings fail the most (for debugging).
   */
  static vector<int> ofail, o1fail;

  void fail(int i ) const;

  /**
   * Exception class to signal bad kinematics.
   */
  struct InteractionKinematicException: public Exception {};

};

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
  inline DipoleXSec()
    : theRMax(0.0*InvGeV), theInteraction(0), sinFunction(0), usePartonicInteraction(false),
      theIntOrdering(0), theRecoilReduction(0), checkOffShell(true) {}

  /**
   * The copy constructor.
   */
  inline DipoleXSec(const DipoleXSec & x)
    : HandlerBase(x), theRMax(x.theRMax), theInteraction(x.theInteraction),
      sinFunction(x.sinFunction), usePartonicInteraction(x.usePartonicInteraction),
      theIntOrdering(x.theIntOrdering),theRecoilReduction(x.theRecoilReduction),
      checkOffShell(x.checkOffShell) {}

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
   * Calculate the scattering probability for the two given dipoles
   * using the given ImpactParameters.
   */
  virtual DipoleInteraction fij(const ImpactParameters &  b, Dipole &, Dipole &,
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
   * Return false if an interaction would be kinematically forbidden.
   */
  virtual bool kinematicsVeto(const DipoleInteraction &) const;

  /**
   * Calculate the total scattering probability for the two given
   * dipole systems using the given ImpactParameters.
   */
  virtual double sumf(const DipoleState &, const DipoleState &,
		      const ImpactParameters &) const;

  /**
   * Calculate the total scattering probability for the two given
   * dipole systems using the given ImpactParameters.
   */
  virtual double sumf(const ImpactParameters &,
		      const DipoleState &, const DipoleState &) const;

  /**
   * Calculate the total scattering probability all dipole pairs 
   * the two given dipole systems using the given ImpactParameters.
   */
  virtual FList flist(const DipoleState &, const DipoleState &,
		      const ImpactParameters &) const;

  /**
   * Calculate the total scattering probability all dipole pairs 
   * the two given dipole systems using the given ImpactParameters.
   */
  virtual DipoleInteraction::List flist(const ImpactParameters &,
				const DipoleState &, const DipoleState &) const;


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
   * Returns the recoils the interaction would like to give the 4 involved partons.
   * Assumes the states are not yet rotated.
   */
  virtual InteractionRecoil
  recoil(const DipoleInteraction &) const;

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
   * creates and initialises a RealInteraction with realpartons, ranges and max.
   */
  virtual RealInteraction
  initialiseInteraction(const DipoleInteraction & di,
			RealPartonStatePtr lrs, RealPartonStatePtr rrs) const;

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
   * Does the transverse recoil on the real partons with the pT supplied in recs.
   * Does the p+ alt. p- transfer of neededPlus nad neededMinus to put the
   * other state on shell. If not consistent, return false.
   * Also updates p_mu for the 4 effective partons. Assumes all are rightmoving.
   */
  virtual bool doInteraction(InteractionRecoil recs, 
			     const DipoleInteraction::List::const_iterator inter,
			     RealPartonStatePtr lrs, RealPartonStatePtr rrs,
			     pair<pair<bool, bool>, pair<bool, bool> > doesInt) const;


  /**
   * reconnects the colour flow in an interaction, taking care that rescatterings are
   * treated properly.
   */
  bool reconnect(const DipoleInteraction &di) const {
    return reconnect(di.dips.first, di.dips.second);
  }

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
   * Mark two dipoles as having interacted with eachother.
   */
  virtual void interact(Dipole & d1, Dipole & d2) const;

  //@}

  /**
   * Return the value of the sine-interaction strength.
   */
  double fSinFn(const Parton::Point & rho1, const Parton::Point & rho2,
		const TransverseMomentum & pt) const;

  /**
   * Return the interaction
   **/
  inline int interaction() const {
    return theInteraction;
  }

  /**
   * Flag determining if only one parton in each dipole is considered
   * interacting or both.
   */
  inline bool partonicInteraction() const {
    return usePartonicInteraction;
  }

  /**
   * The confinement scale.
   */
  InvEnergy rMax() const;

  /**
   * for debugging
   */
  mutable int Nfij;
  mutable int nIAccepted;
  mutable int nIBelowCut;
  mutable int nIPropFail;
  mutable int nIKineFail;
  mutable int nIOrdering;
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
  inline virtual IBPtr clone() const {
    return new_ptr(*this);
  }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {
    return new_ptr(*this);
  }
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
   * Flag determining which approximation to the sine-functions in the
   * interaction strength to use.
   */
  int sinFunction;

  /**
   * Flag determining if only one parton in each dipole is considered
   * interacting or both.
   */
  bool usePartonicInteraction;

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

#endif /* DIPSY_DipoleXSec_H */
