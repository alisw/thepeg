// -*- C++ -*-
#ifndef DIPSY_DipoleAbsorber_H
#define DIPSY_DipoleAbsorber_H
//
// This is the declaration of the DipoleAbsorber class.
//

#include "DipoleAbsorber.fh"
#include "ThePEG/Handlers/HandlerBase.h"
#include "Ariadne/DIPSY/DipoleState.fh"
#include "Ariadne/DIPSY/Parton.fh"
#include "Ariadne/DIPSY/Dipole.fh"

namespace DIPSY {

using namespace ThePEG;

/**
 * DipoleAbsorber is a base class to be used for models describing how
 * non-interacted dipoles will be reabsorbed. The only function to be
 * overridden by sub-classes is called reabsorb.
 *
 * @see \ref DipoleAbsorberInterfaces "The interfaces"
 * defined for DipoleAbsorber.
 */
class DipoleAbsorber: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleAbsorber();

  /**
   * The destructor.
   */
  virtual ~DipoleAbsorber();
  //@}

public:

  /** @name Virtual functions to be overridden in sub-classes. */
  //@{
  /**
   * Reabsorb dipoles in the given state.
   */
  virtual void reabsorb( DipoleState &) const = 0;

  /**
   * Swing so that all non participating nucleons are just the original triangle.
   */
  virtual void isolateNonParticipating ( DipoleState & ) const;

  /**
   * Adds recoils from dipoles that has emitted, trying to account
   * for emission history.
   */
  virtual void hiddenRecoil(DipoleState & ) const = 0;

  /**
   * Absorb the provided parton, sharing p_mu according to distance.
   * swing determines if a swing should be forced if the parton is in a loop.
   */
  virtual void absorbParton ( tPartonPtr p, DipoleState & ) const = 0;

  /**
   * Removes the parton and recouples the colow flow,
   * without transfering its momentum to other partons.
   * Does not swing loops, that has to be handled before calling this.
   */
  virtual void removeParton ( tPartonPtr p ) const;

  /**
   * Absorb the virtual partons that has no interacting children.
   */
  virtual void absorbVirtualPartons ( DipoleState & ) const;

  /**
   * Swings the loop that the dipole /a d is part of.
   */
  virtual void swingLoop( DipolePtr d, DipoleState &) const;

  /**
   * Swings the two dipoles.
   */
  virtual bool swing ( DipolePtr, DipolePtr ) const;

  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleAbsorber & operator=(const DipoleAbsorber &);

};

}

#endif /* DIPSY_DipoleAbsorber_H */
