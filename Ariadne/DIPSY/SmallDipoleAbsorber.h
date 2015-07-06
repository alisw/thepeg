// -*- C++ -*-
#ifndef DIPSY_SmallDipoleAbsorber_H
#define DIPSY_SmallDipoleAbsorber_H
//
// This is the declaration of the SmallDipoleAbsorber class.
//

#include "Ariadne/DIPSY/DipoleAbsorber.h"
#include "DipoleState.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * The SmallDipoleAbsorber class implements a model describing how
 * non-interacted dipoles will be reabsorbed. The model is based on
 * the removal of small dipoles.
 *
 * @see \ref SmallDipoleAbsorberInterfaces "The interfaces"
 * defined for SmallDipoleAbsorber.
 */
class SmallDipoleAbsorber: public DipoleAbsorber {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SmallDipoleAbsorber();

  /**
   * The destructor.
   */
  virtual ~SmallDipoleAbsorber();
  //@}

public:

  /** @name Virtual functions to be overridden. */
  //@{
  /**
   * Reabsorb dipoles in the given state.
   */
  virtual void reabsorb( DipoleState &) const;

  // /**
  //  * Swing so that all non participating nucleons are just the original triangle.
  //  */
  // virtual void isolateNonParticipating ( DipoleState & ) const;
  //@}

public:

  /**
   * Adds recoils from dipoles that has emitted, trying to account
   * for emission history.
   */
  void hiddenRecoil(DipoleState &) const;

  /**
   * Absorb the provided parton, sharing p_mu according to distance.
   * swing determines if a swing should be forced if the parton is in a loop.
   */
  virtual void absorbParton (tPartonPtr p, DipoleState & ) const;

  /**
   * Swings the two dipoles.
   */
  virtual bool swing (tDipolePtr, tDipolePtr ) const;

  /**
   * Absorb the nonordered partons according to sorted.
   */
template <typename Sort>
  void absorbNonOrderedPartons ( DipoleState &, Sort sort ) const;

  /**
   * Sets the parton as sorted, and recursively also all its neighbours.
   */
template < typename Sort >
  void setOrdered ( tPartonPtr p, Sort sort ) const;

  /**
   * Confirm the dipole as real, and recursively also all its larger neighbours.
   */
  void setDGLAPsafe ( tDipolePtr dip ) const;

  /**
   * Decides if a dipole of length r1 can be resolved by a dipole of length r2.
   */
  bool checkDGLAPsafe ( InvEnergy r1, InvEnergy r2 ) const;

  /**
   * compoensate for the 1/r^2 vs 1/r^4 weights by absorbing small dipoles created
   * with too high weight.
   */
  void absorbSmallDipoles ( DipoleState & ) const;

  /**
   * Set the parton as y ordered, and recursively also all its ordered neighbours.
   */
  void setYOrdered ( tPartonPtr p ) const;

  /**
   * Set the parton as p- ordered, and recursively also all its ordered neighbours.
   */
  void setMinusOrdered ( tPartonPtr p ) const;

  /**
   * Absorb the parton, and all partons in the direction of smaller dipoles.
   * Starts with the smallest dipole at the end of the chain.
   */
  void DGLAPabsorb ( tPartonPtr p, DipoleState & ) const;

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

private:

  /**
   * Recoils the two partons, as from an emission. pT and y changed.
   * p+/p- conserved for right/left moving particles respectively. y changes.
   * ymax is evolved rapidity of the state, to know when the partons recoil
   * too far in rapidity.
   */
  void recoil(tPartonPtr p1, tPartonPtr p2, double ymax) const;

  /**
   * Marks the valence partons as interacted, that is non virtual.
   */
  void realValence(DipoleState &) const;

  /**
   * Finds the largest non DGLAPsafe dipole among /a dipoles that is
   * not in a loop that should be absorbed. Returns an empty pointer
   * if no nonsafe dipoles that shouldnt be absorb were found.
   */
  DipolePtr largestNonsafeDipole(list<DipolePtr> dipoles, DipoleState &) const;

  /**
   * Absorb the loop that the dipole /a d is part of.
   */
  void absorbLoop( DipolePtr d, DipoleState &) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SmallDipoleAbsorber & operator=(const SmallDipoleAbsorber &);

};

}

#endif /* DIPSY_SmallDipoleAbsorber_H */
