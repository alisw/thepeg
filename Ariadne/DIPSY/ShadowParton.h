// -*- C++ -*-
#ifndef DIPSY_ShadowParton_H
#define DIPSY_ShadowParton_H
//
// This is the declaration of the ShadowParton class.
//

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Vectors/Transverse.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/UseRandom.h"
#include "Ariadne/Config/CloneBase.h"
#include "ShadowParton.fh"
#include "Parton.h"
#include "Dipole.fh"
#include "DipoleXSec.fh"
#include "Sum20Momentum.h"

namespace DIPSY {

using namespace ThePEG;

class ImpactParameters;

/**
 * Here is the documentation of the ShadowParton class.
 */
class ShadowParton: public Ariadne5::CloneBase {

public:

  /**
   * Typedef for position in transverse coordinate space.
   */
  typedef Transverse<InvEnergy> Point;

  /**
   * A pair of partons.
   */
  typedef pair<tSPartonPtr,tSPartonPtr> tSPartonPair;

  /**
   * Helper class for traversing along a propagator to an interaction
   */
  struct Propagator {

    /**
     * Constructor.
     */
    Propagator(): colpos(Constants::MaxEnergy), colneg(ZERO), colptp(ZERO),
		  acopos(Constants::MaxEnergy), aconeg(ZERO), acoptp(ZERO),
		  fail(false) {}

    /**
     * Copy constructor.
     */
    Propagator(const Propagator & k)
      : p(k.p), colpos(k.colpos), colneg(k.colneg), colptp(k.colptp),
	acopos(k.acopos), aconeg(k.aconeg), acoptp(k.acoptp), fail(k.fail) {}

    /**
     * Construct from momentum.
    p */
    Propagator(const LorentzMomentum & mom)
      : p(mom), colpos(Constants::MaxEnergy), colneg(ZERO), colptp(ZERO),
	acopos(Constants::MaxEnergy), aconeg(ZERO), acoptp(ZERO), fail(false) {}

    /**
     * The Momentum of the propagator.
     */
    LorentzMomentum p;

    /**
     * The positive light-cone momentum of the previous emission on
     * the colour line.
     */
    Energy colpos;

    /**
     * The negative light-cone momentum of the previous emission on
     * the colour line.
     */
    Energy colneg;

    /**
     * The pt of the propagator before the previous emission on the
     * colour line (ie. previous emission with an inti-colour sibling).
     */
    Energy2 colptp;

    /**
     * The positive light-cone momentum of the previous emission on
     * the anti-colour line (ie. previous emission with a colour sibling).
     */
    Energy acopos;

    /**
     * Thenegative light-cone momentum of the previous emission on the
     * anti-colour line.
     */
    Energy aconeg;

    /**
     * The pt of the propagator before the previous emission on the
     * anti-colour line.
     */
    Energy2 acoptp;

    /**
     * Indicate that the propagator has failed.
     */
    bool fail;

  };

  /**
   * Helper class to temporarily protect a propagator.
   */
  struct Lock {

    /**
     * The only constructor. Locking the propagator to a given parton.
     */
    Lock(tPartonPtr p): sp(p->shadow()) {
      sp->lock();
    }

    /**
     * The destructor, unlocking the propagator.
     */
    ~Lock() {
      sp->unlock();
    }

    /**
     * The parton to which the propagator should be locked.
     */
    tSPartonPtr sp;

  };


public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ShadowParton()
    : thePlus(ZERO), theMass(-1*GeV), theMinus(0*GeV),
      theY(0), theY0(0), theFlavour(ParticleID::g),
      theEmissionFactor(ZERO), theRes(ZERO), isLocked(false),
      hasInteracted(false), isOnShell(false), isValence(false),
      hasColourSibling(false), memememe(false) {}


  /**
   * Construct from ordinary parton.
   */
  ShadowParton(Parton & p);

  /**
   * The destructor.
   */
  virtual ~ShadowParton();

  /**
   * Create a valence shadow parton from a given ordinary parton.
   */
  static SPartonPtr createValence(Parton & p);

  /**
   * Easy access to alphaS.
   */
  static double alphaS(InvEnergy2 r2);

  /**
   * Setup shadow partons in an emission.
   */
  void setupEmission(Parton & emitter, Parton & produced, Parton & recoiler);

  /**
   * Setup pointers to and from neighboring partons.
   */
  void setupParton();

  //@}

protected:

  /** @name The virtual functions to be overridden in sub-classes. */
  //@{
  /**
   * Return a simple clone of this object. Should be implemented as
   * <code>return new_ptr(*this);</code> by a derived class.
   */
  virtual Ariadne5::ClonePtr clone() const;

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

  /**
   * Return the DipoleState to which this parton belongs. If it was
   * originally belonged to one state which was subsequently merged
   * into another, the parton will belong to the latter state.
   */
  DipoleState & dipoleState() const;

  /**
   * finds the pTscale through the dipolestates emitter.
   */
  double pTScale() const;

  /**
   * Calculate the squared transverse distance to the given parton.
   */
  inline InvEnergy2 dist2(const ShadowParton & p) const {
    return (position() - p.position()).pt2();
  }

  /**
   * Calculate the transverse distance to the given parton.
   */
  inline InvEnergy dist(const ShadowParton & p) const {
    return sqrt(dist2(p));
  }

  /**
   * Produce a ThePEG::Particle corresponding to this parton.
   */
  PPtr produceParticle() const;

  /**
   * The mass of this parton.
   */
  Energy mass() const;

  /**
   * The transverse momentum squared of this parton.
   */
  Energy2 pt2() const {
    return pT().pt2();
  }

  /**
   * The transverse mass squared of this parton.
   */
  Energy2 mt2() const {
    return pt2() + sqr(mass());
  }

  /**
   * The transverse mass of this parton in the emission.
   */
  Energy2 mt02() const {
    return pT0().pt2() + sqr(mass());
  }

  /**
   * The transverse mass of this parton in the emission.
   */
  Energy mt0() const {
    return sqrt(mt02());
  }

  /**
   * The transverse momentum of this parton.
   */
  Energy pt() const {
    return sqrt(pt2());
  }

  /**
   * The transverse momentum of this parton.
   */
  Energy mt() const {
    return sqrt(mt2());
  }

  /**
   * The final-state momentum of this particle.
   */
  LorentzMomentum momentum() const {
    return lightCone(plus(), mt2()/plus(), pT());
  }

  /** @name Simple access functions. */
  //@{
  /**
   * Get the position in impact parameter space.
   */
  inline const Point & position() const {
    return original()->position();
  }

  /**
   * Get the positive light-cone momentum.
   */
  inline Energy plus() const {
    return thePlus;
  }

  /**
   * Get the original positive light-cone momentum.
   */
  inline Energy plus0() const {
    return pT0().pt()*exp(-y0());
  }

  /**
   * Get the transverse momentum.
   */
  inline TransverseMomentum pT() const {
    return thePT;
  }

  /**
   * Get the transverse momentum.
   */
  inline TransverseMomentum pT0() const {
    return thePT0;
  }

  /**
   * Get the accumulated negative light-cone momentum deficit.
   */
  inline Energy minus() const {
    return theMinus;
  }

  /**
   * Get the rapidity.
   */
  inline double y() const {
    return theY;
  }

  /**
   * Get the rapidity generated in the emission.
   */
  inline double y0() const {
    return theY0;
  }

  /**
   * Get the flavour of this parton.
   */
  inline long flavour() const {
    return theFlavour;
  }

  /**
   * Get the original DIPSY parton.
   */
  tPartonPtr original() const {
    return theOriginal;
  }

  /**
   * Get the previous instance of this parton befor it was emitted or
   * swinged.
   */
  inline tSPartonPtr previous() const {
    return thePrevious;
  }

  /**
   * Get the parent parton if this was emitted.
   */
  inline tSPartonPtr parent() const {
    return theParent;
  }

  /**
   * Get the sibling of this parton if any.
   */
  inline tSPartonPtr sibling() const {
    return theSibling;
  }
		      
  /**
   * Get the next version of this parton if emitted or swinged.
   */
  inline tSPartonPtr next() const {
    return theNext;
  }

  /**
   * Get the first version of this parton if emitted or swinged.
   */
  tSPartonPtr initial();
		      
  /**
   * Get the first (grand)mother of this parton.
   */
  tSPartonPtr valenceMother();
		      
  /**
   * Get the last version of this parton if emitted or swinged.
   */
  tSPartonPtr last();
		      
  /**
   * Get the last version of this parton if emitted or swinged.
   */
  tcSPartonPtr last() const;
		      
  /**
   * Get the parton this has emitted if any.
   */
  inline tSPartonPtr child() const {
    return theChild;
  }

  /**
   * Get the mother parton, either previous() or parent().
   */
  inline tSPartonPtr mother() const {
    return parent()? parent(): previous();
  }

  /**
   * Indicate that the corresponding emission is on a protected
   * propagator. Any other propagator must assume that this emission
   * must be resolved.
   */
  inline bool locked() const {
    return isLocked;
  }

  /**
   * Indicate that the propagator leading up to this parton is protected.
   */
  void lock();

  /**
   * Indicate that the propagator leading up to thi parton is no
   * longer protected.
   */
  void unlock();

  /**
   * Indicate if this parton has interacted.
   */
  inline bool interacted() const {
    return hasInteracted;
  }

  /**
   * Return true if the sibling is on the colour side of this parton.
   */
  inline bool colourSibling() const {
    return hasColourSibling;
  }

  /**
   * Return if this shadow was the root of an interacting propagator.
   */
  inline bool interactionRoot() const {
    return theInteractionRoots.size() > 0;
  }

  /**
   * Return if this shadow was the root of an interacting propagator.
   */
  inline bool hasInteractionRoot(int i) const {
    return theInteractionRoots.find(i) != theInteractionRoots.end();
  }

  /**
   * Return if this shadow was directly involved in the interaction.
   */
  inline int interacting() const {
    return theInteractions.size() > 0;
  }

  /**
   * Return if this parton is on shell.
   */
  inline bool onShell() const {
    return isOnShell;
  }

  /**
   * Set if this parton is on shell.
   */
  inline void onShell(bool b) {
    isOnShell = b;
  }

  /**
   * Indicate that this parton is the root of an interaction.
   */
  inline void interactionRoot(int i) {
    theInteractionRoots.insert(i);
  }

  /**
   * Set if this shadow is directly involved in an interaction.
   */
  inline void interacting(int i) {
    theInteractions.insert(i);
  }

  /**
   * Return if this parton is a valence parton.
   */
  inline bool valence() const {
    return isValence;
  }

  /**
   * Set if this parton is a valence parton.
   */
  inline void valence(bool b) {
    isValence = b;
  }

  /**
   * Set the positive light-cone momentum.
   */
  inline void plus(Energy x) {
    thePlus = x;
  }

  /**
   * Set the transverse momentum.
   */
   inline void pT(const TransverseMomentum & x) {
     thePT = x;
   }

  /**
   * Set the momentum using transverse momentum and positive
   * light-cone momentum. Set also other components assuming real.
   */
  void pTplus(const TransverseMomentum & qt, Energy qp);

  /**
   * Reset the momentum to its original value in the evolution.
   */
  void resetMomentum0() {
    pTplus(pT0(), mt0()/exp(y0()));
  }

  /**
   * Set the transverse momentum.
   */
   inline void pT0(const TransverseMomentum & x) {
     thePT0 = x;
   }

  /**
   * Set the accumulated negative light-cone momentum deficit.
   */
  inline void minus(Energy x) {
    theMinus = x;
  }

  /**
   * Set the rapidity.
   */
  inline void y(double x) {
    theY = x;
  }

  /**
   * Set the rapidity.
   */
  inline void y0(double x) {
    theY0 = x;
  }

  /**
   * Set the flavour of this parton.
   */
  inline void flavour(long x) {
    theFlavour = (x == 0? ParticleID::g: x);
  }

  /**
   * Set the parent parton.
   */
  inline void parent(tSPartonPtr p) {
    theParent = p;
  }

  /**
   * Indicate that this parton has interacted and mark all parent and
   * original partons if needed.
   */
  void interact();

  /**
   * The emission factor for this parton if emitted or emitter.
   */
  InvEnergy2 emissionFactor() const {
    return theEmissionFactor;
  }

  /**
   * The resolusion scale factor for this parton if emitted or emitter.
   */
  InvEnergy2 res() const {
    return theRes;
  }

  /**
   * Return true if this parton and its sibling are resolved.
   */
  bool resolved(InvEnergy2 r2) const {
    return ( ( sibling() && ( sibling()->locked()
			      || emissionFactor() > alphaS(r2)*r2 ) ) );
  }

  /**
   * Return true if the sibling must be put on-shell even if it was not
   * resolved, because it will partake in a subsequent interaction or is
   * a valence.
   */
  bool forceEmission() const {
    return sibling() && ( sibling()->valence()
			  || sibling()->interactionRoot() );
  }

  /**
   * Go back in the history of this shadow parton and find the
   * original parton which would be resolved at a given size or would
   * be the interaction root. Always resolve emissions involving \a
   * stopp.
   */
  tSPartonPtr resolveInteraction(InvEnergy2 r2, tPartonPtr stopp) {
    Lock safeguard(stopp);
    return resolveInteraction(r2);
  }

  /**
   * Go back in the history of this shadow parton and find the
   * original parton which would be resolved at a given size or would
   * be the interaction root.
   */
  tSPartonPtr resolveInteraction(InvEnergy2 r2);

  /**
   * Go back in the history of this shadow parton and find the
   * original parton which would be resolved at a given size. Always
   * resolve emissions involving \a stopp.
   */
  tSPartonPtr resolve(InvEnergy2 r2, tPartonPtr stopp) {
    Lock safeguard(stopp);
    return resolve(r2);
  }

  /**
   * Go back in the history of this shadow parton and find the
   * original parton which would be resolved at a given size. Always
   * resolve emissions involving \a stopp. (const version).
   */
  //  tcSPartonPtr resolve(InvEnergy2 r2, tPartonPtr stopp) const;

  /**
   * Go back in the history and flag all partons which are resoved by
   * the given size to be set on-shell. Always
   * resolve emissions involving \a stopp.
   */
  //  void flagOnShell(InvEnergy2 r2, tPartonPtr stopp);

  /**
   * If this parton is not on shell, find an emitted parton
   * with the same colour lines (modulo swing) that is.
   */
  tSPartonPtr findFirstOnShell();

  /**
   * If this parton is not on shell, find an emitted parton
   * with the same colour lines (modulo swing) that is.
   */
  tSPartonPtr findSecondOnShell();

  /**
   * Flag on-shell if \a mode >= 0 and copy to original parton if \a
   * mode > 0.
   */
  void setOnShell(int mode);

  /**
   * If prevoíously set on-shell, clear the flag.
   */
  void unsetOnShell();

  /**
   * Go forward in the evolution and reset all interaction flags.
   */
  void reset();

  /**
   * Go forward in the evolution and reset all interacted flags.
   */
  void resetInteracted();

  /**
   * Set the momenta in an emission, assuming the incoming momentum \a p.
   */
  bool setEmissionMomenta(const LorentzMomentum & p, bool valencefix);

protected:

  /**
   * Go back in the history of this shadow parton and find the
   * original parton which would be resolved at a given size. Always
   * resolve emissions involving \a stopp.
   */
  tSPartonPtr resolve(InvEnergy2 r2);

public:

  /**
   * Go back in the history and find the momentum of this incoming
   * parton. Always resolve emissions involving \a stopp. If \a mode
   * >= 0 flag partons which would be put on shell. if \a mode > 0
   * also propagate the momentum to the original() parton.
   */
  Propagator propagator(InvEnergy2 r2, tPartonPtr stopp, int mode) {
    Lock safeguard(stopp);
    return propagator(r2, mode);
  }
  Propagator propagator(InvEnergy2 r2, int mode);

  Propagator setMomentum(const Propagator & prop) {
    plus(prop.p.plus());
    pT(TransverseMomentum(prop.p.x(), prop.p.y()));
    return prop;       
  }

  void checkMomentum(Sum20Momentum & sum20,
		     const ImpactParameters * b = 0) const;

  /**
   * Indicate that this parton should be tested for an
   * interaction. Insert the corresponding interaction object at the
   * correct place in the chain. The primary interaction is always at
   * one of the incoming shadows. Secondary interactions are inserted
   * at a parton which has been or should be made real by a previous
   * emission. If this parton already had an interaction, simply leave
   * it as it is.
   */
  void insertInteraction(int i);

  /**
   * Indicate that a previously tested interaction is accepted. This
   * is done by marking all propagators interacted() in the chain down
   * to the one where the interaction was found.
   */
  void acceptInteraction(int i);

  /**
   * Make this and all its anceseters incoming.
   */
  void makeIncoming();

  /**
   * Indicate that the interaction that has been tested was rejected
   * by simply removing the interaction root inserted by
   * insertInteraction().
   */
  void rejectInteraction(int i);
  

  /**
   * Change direction before producing collision.
   */
  void mirror(double yframe);

  /**
   * Translate according to the impagt parameter.
   */ 
  void translate(const ImpactParameters & b);

  /**
   * Check if the generated emission fails the ordering requirement.
   */
  bool orderfail(const Propagator & prop) const;

  /**
   * Setup the continuing propagater with the generated emission.
   */
  Propagator setup(Propagator & prop) const;

  /**
   * Prints all the member variables to cerr.
   */
  void debugme();

  /**
   * Debug the structure of the shadow tree.
   */
  void debugTree(string indent);

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
   * The original DIPSY parton
   */
  tPartonPtr theOriginal;

  /**
   * The positive light-cone momentum.
   */
  Energy thePlus;

  /**
   * The transverse momentum.
   */
  TransverseMomentum thePT;

  /**
   * The transverse momentum generated in the emission.
   */
  TransverseMomentum thePT0;

  /**
   * The mass of this parton.
   */
  mutable Energy theMass;

  /**
   * The accumulated negative light-cone momentum deficit.
   */
  Energy theMinus;

  /**
   * The rapidity.
   */
  double theY;

  /**
   * The rapidity generated in the emission.
   */
  double theY0;

  /**
   * The flavour of this parton.
   */
  long theFlavour;

  /**
   * The state of this parton, if any, before the swing or the emission of this
   * parton's sibling. Is aways null if theParent exists.
   */
  tSPartonPtr thePrevious;

  /**
   * The parent parton if this has emitted or swinged. Is aways null if
   * thePrevious exists.
   */
  tSPartonPtr theParent;

  /**
   * The sibling of this parton if any.
   */
  tSPartonPtr theSibling;

  /**
   * The state of this parton after the emission of the child. Is null
   * if this parton has not emitted.
   */
  SPartonPtr theNext;

  /**
   * The child this parton has emitted, if any.
   */
  SPartonPtr theChild;

  /**
   * The emission factor for this parton if emitted or emitter.
   */
  InvEnergy2 theEmissionFactor;

  /**
   * The resolution factor for this parton if emitted or emitter.
   */
  InvEnergy2 theRes;

  /**
   * Indicate that the corresponding emission is on a protected
   * propagator. Any other propagator must assume that this emission
   * must be resolved.
   */
  bool isLocked;

  /**
   * Indicate if this parton has interacted.
   */
  bool hasInteracted;

  /**
   * Indicate if this parton is on shell, that is if it has
   * been supplied with the neccesary p+ or p- for left or rightmoving
   * particles respectively.
   */
  bool isOnShell;

  /**
   * This is root of an interaction chain.
   */
  set<int> theInteractionRoots;

  /**
   * The interactions this was directly involved with.
   */
  set<int> theInteractions;

  /**
   * Indicate if this parton is a valence parton.
   */
  bool isValence;

  /**
   * Indicate if this partons sibling is on the colour side (true) or
   * the anti-colour side (false).
   */
  bool hasColourSibling;

protected:

  /**
   * Exception class for bad kinematics.
   */
  struct ShadowPartonKinematicsException: public Exception {};

public:

  bool memememe;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShadowParton & operator=(const ShadowParton &);

};

}

#endif /* DIPSY_ShadowParton_H */
