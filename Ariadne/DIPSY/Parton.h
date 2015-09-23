// -*- C++ -*-
#ifndef DIPSY_Parton_H
#define DIPSY_Parton_H
//
// This is the declaration of the Parton class.
//

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Vectors/Transverse.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Ariadne/Config/CloneBase.h"
#include "Parton.fh"
#include "ShadowParton.fh"
#include "Dipole.fh"

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the Parton class.
 */
class Parton: public Ariadne5::CloneBase {

public:

  friend class RealParton;

  /**
   * Typedef for position in transverse coordinate space.
   */
  typedef Transverse<InvEnergy> Point;

  /**
   * A pair of partons.
   */
  typedef pair<tPartonPtr,tPartonPtr> tPartonPair;

  /**
   * A pair of dipoles.
   */
  typedef pair<tDipolePtr,tDipolePtr> tDipolePair;

  /**
   * The Dipole is a friend.
   */
  friend class Dipole;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Parton()
    : thePosition(Point()), thePlus(0.0*GeV), theOriginalY(0.0),
      thePT(TransverseMomentum()),
      theValencePT(TransverseMomentum()),theValencePlus(ZERO),
      theMinus(0*GeV), theY(0), isRightMoving(true), theFlavour(ParticleID::g),
      theParents(tPartonPair()), theDipoles(tDipolePair()), hasInteracted(false),
      isOnShell(false), isValence(false), theMass(-1*GeV), theShadow(tSPartonPtr()) {}

  /**
   * The copy constructor.
   */
  inline Parton(const Parton & x)
    : thePosition(x.thePosition), thePlus(x.thePlus),theOriginalY(x.theOriginalY), 
      thePT(x.thePT), theValencePT(x.theValencePT),theValencePlus(x.theValencePlus),
      theMinus(x.theMinus), theY(x.theY), isRightMoving(x.isRightMoving),
      theFlavour(x.theFlavour),theParents(x.theParents),
      theChildren(x.theChildren), theDipoles(x.theDipoles),
      hasInteracted(x.hasInteracted), isOrdered(x.isOrdered), isOnShell(x.isOnShell),
      isValence(x.isValence), theMass(x.theMass), theShadow(x.theShadow) {}

  /**
   * The destructor.
   */
  virtual ~Parton();
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
  inline InvEnergy2 dist2(const Parton & p) const {
    return (position() - p.position()).pt2();
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
   * The scalar transverse momentum squared of this parton.
   */
  Energy2 pt2() const {
    return pT().pt2();
  }

  /**
   * The scalar transverse mass squared of this parton.
   */
  Energy2 mt2() const {
    return pT().pt2() + sqr(mass());
  }

  /**
   * The scalar transverse momentum of this parton.
   */
  Energy pt() const {
    return sqrt(pt2());
  }

  /**
   * The scalar transverse mass of this parton.
   */
  Energy mt() const {
    return sqrt(mt2());
  }

  /**
   * The final-state momentum of this particle.
   */
  LorentzMomentum momentum() const {
    return lightCone(plus(), minus(), pT());
  }

  /** @name Simple access functions. */
  //@{
  /**
   * Get the position in impact parameter space.
   */
  inline const Point & position() const {
    return thePosition;
  }

  /**
   * Get the positive light-cone momentum.
   */
  inline Energy plus() const {
    return thePlus;
  }

  /**
   * Get the transverse momentum.
   */
  inline TransverseMomentum pT() const {
    return thePT;
  }

  /**
   * Get the transverse valence momentum. 0 if not a valence parton.
   */
  inline TransverseMomentum valencePT() const {
    if (isValence)
      return theValencePT;
    return TransverseMomentum();
  }

  /**
   * Get the valence energy. 0 if not a valence parton.
   */
  inline Energy valencePlus() const {
    if (isValence)
      return theValencePlus;
    return ZERO;
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
   * Get the flavour of this parton.
   */
  inline long flavour() const {
    return theFlavour;
  }

  /**
   * Get the parent partons.
   */
  inline tPartonPair parents() const {
    return theParents;
  }

  /**
   * Get the main parent.
   */
  inline tPartonPtr mainParent() const {
    return theMainParent;
  }

  /**
   * Get the children partons.
   */
  inline const set<tPartonPtr> & children() const {
    return theChildren;
  };

  /**
   * Remove a child
   **/
  inline void removeChild(PartonPtr child) {
    theChildren.erase(child);
  };

  /**
   * Remove all children
   **/
  inline void removeAllChildren() {
    theChildren.clear();
  };

  /**
   * Get the connecting dipoles.
   */
  inline tDipolePair dipoles() const {
    return theDipoles;
  }

  /**
   * Indicate if this parton has interacted.
   */
  inline bool interacted() const {
    return hasInteracted;
  }

  /**
   * Indicate if this parton is ordered.
   */
  inline bool ordered() const {
    return isOrdered;
  }

  /**
   * Set if this parton is ordered.
   */
  inline void ordered(bool b) {
    isOrdered = b;
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
   * Return if this parton is a valence parton.
   */
  inline bool valence() const {
    return isValence;
  }

  /**
   * Returns true if the parton is absorbed.
   */
  bool absorbed() const;

  /**
   * Set if this parton is a valence parton.
   */
  inline void valence(bool b) {
    isValence = b;
  }

  /**
   * Set the position in impact parameter space.
   */
  inline void position(const Point & x) {
    thePosition = x;
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
   inline void pT(TransverseMomentum x) {
     thePT = x;
   }

  /**
   * Set the transverse valence momentum, and set isValence to true.
   */
   inline void valencePT(TransverseMomentum x) {
     isValence = true;
     theValencePT = x;
   }

  /**
   * Set the valence energy.
   */
   inline void valencePlus(Energy E) {
     theValencePlus = E;
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
   * Set the flavour of this parton.
   */
  inline void flavour(long x) {
    theFlavour = (x == 0? ParticleID::g: x);
  }

  /**
   * Set the parent partons.
   */
  inline void parents(tPartonPair x) {
    theParents = x;
  }

  /**
   * Set the main parent parton.
   */
  inline void mainParent(tPartonPtr p) {
    theMainParent = p;
  }

  /**
   * Set the connecting dipoles.
   */
  inline void dipoles(tDipolePair x) {
    theDipoles = x;
  }

  /**
   * Get the original emission rapidity.
   */
  inline const double oY() const {
    return theOriginalY;
  }

  /**
   * Set the original emission rapidity.
   */
  inline void oY(double y) {
    theOriginalY = y;
  }

  /**
   * Indicate that this parton has interacted. If \a onegluon is true,
   * only the main parent is recursively flagged as interacted.
   */
  void interact(bool onegluon = false);

  /**
   * Returns true if the parton comes from a rightmoving state.
   */
  inline bool rightMoving() const {
    return isRightMoving;
  }

  /**
   * Set the direction of the original state.
   */
  inline void rightMoving(bool rM) {
    isRightMoving = rM;
  }

  /**
   * Get and set the number of the parton.
   */
  inline int number() {
    return theNumber;
  }

  inline void number(int N) {
    theNumber = N;
  }

  /**
   * Checks if the parton is colour connected (in many steps if 
   * needed) to an interacting dipole.
   */
  bool inInteractingLoop() const;

  /**
   * counts the number of on shell partons in this partons colour chain.
   */
  int nOnShellInChain() const;

  /**
   * Checks if the parton is part of on of its states original dipoles.
   */
  bool valenceParton() const;

  /**
   * Checks if the parton is in a colour chain containing a valenceparton.
   */
  bool inValenceChain() const;

  /**
   * Returns true if the dipole the parton was emitted from is the result of a swing.
   **/
  bool swingedEmission() const;

  /**
   * Updates pT according to the colour neighbours, and changes y and
   * plus/minus for right/left moving. minus/plus unchanged.
   */
  void updateMomentum();

  /**
   * Updates y and p- from plus and pT.
   */
  inline void updateYMinus() {
    theY = log(thePT.pt()/thePlus);
    theMinus = thePT.pt()*exp(theY);
  }

  inline void mirror(double y0) {
    y(2.0*y0 - y());
    swap(thePlus, theMinus);
    rightMoving(!rightMoving());
  }

  /**
   * Prints all the member variables to cout.
   */
  void coutData();

  /**
   * Returns the recoil this parton would get from the argument parton.
   */
  TransverseMomentum recoil(tPartonPtr) const;

  /**
   * Return the current shadow parton.
   */
  tSPartonPtr shadow() const {
    return theShadow;
  }

  /**
   * Return a new shadow parton.
   */
  void shadow(tSPartonPtr sp) {
    theShadow = sp;
  }

 //@}

public:
  /**
   * debugging
   */
  mutable TransverseMomentum m1pT, m2pT;
  mutable Energy m1plus, m2plus;
  mutable InvEnergy minr1, minr2, maxr1, maxr2;

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
   * The position in impact parameter space.
   */
  Point thePosition;

  /**
   * The positive light-cone momentum.
   */
  Energy thePlus;

  /**
   * The enumeration of this emission.
   */
  double theOriginalY;

  /**
   * The transverse momentum.
   */
  TransverseMomentum thePT;

  /**
   * The transverse momentum and energy a valenceparton got at creation.
   */
  TransverseMomentum theValencePT;
  Energy theValencePlus;

  /**
   * The accumulated negative light-cone momentum deficit.
   */
  Energy theMinus;

  /**
   * The rapidity.
   */
  double theY;

  /**
   * The original direction of the particle.
   */
  bool isRightMoving;

  /**
   * The flavour of this parton.
   */
  long theFlavour;

  /**
   * The parent partons.
   */
  tPartonPair theParents;

  /**
   * The main parent parton.
   */
  tPartonPtr theMainParent;

  /**
   * The children
   **/
  set<tPartonPtr> theChildren;

  /**
   * The connecting dipoles.
   */
  tDipolePair theDipoles;

  /**
   * Indicate if this parton has interacted.
   */
  bool hasInteracted;

  /**
   * Indicate if this parton is ordered with respect to its colour neighbours.
   */
  bool isOrdered;

  /**
   * Indicate if this parton is on shell, that is if it has
   * been supplied with the neccesary p+ or p- for left or rightmoving
   * particles respectively.
   */
  bool isOnShell;

  /**
   * Indicate if this parton is a valence parton.
   */
  bool isValence;

  /**
   * A number that identifies the parton in a dipole state.
   */
  int theNumber;

  /**
   * The mass of this parton.
   */
  mutable Energy theMass;

  /**
   * The current shadow parton.
   */
  tSPartonPtr theShadow;

protected:

  /**
   * Exception class for bad kinematics.
   */
  struct PartonKinematicsException: public Exception {};

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Parton & operator=(const Parton &);

};

}

template <typename T, typename Alloc>
std::vector<T,Alloc> &
expandToFit(std::vector<T,Alloc> & v,
	    typename std::vector<T,Alloc>::size_type indx) {
  if ( indx >= v.size() ) v.resize(indx + 1);
  return v;
}

template <typename T, typename Alloc>
T & forceAt(std::vector<T,Alloc> & v,
	    typename std::vector<T,Alloc>::size_type indx) {
  return expandToFit(v, indx)[indx];
}

#endif /* DIPSY_Parton_H */
