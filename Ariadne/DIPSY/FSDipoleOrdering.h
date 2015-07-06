// -*- C++ -*-
#ifndef DIPSY_FSDipoleOrdering_H
#define DIPSY_FSDipoleOrdering_H
//
// This is the declaration of the FSDipoleOrdering class.
//

#include "Ariadne/DipoleCascade/ReweightBase.h"

namespace DIPSY {

using namespace Ariadne;

/**
 * Here is the documentation of the FSDipoleOrdering class.
 *
 * @see \ref FSDipoleOrderingInterfaces "The interfaces"
 * defined for FSDipoleOrdering.
 */
class FSDipoleOrdering: public Ariadne::ReweightBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FSDipoleOrdering();

  /**
   * The destructor.
   */
  virtual ~FSDipoleOrdering();
  //@}

public:

  /**
   * In addition to the reweight function a final hit/miss veto may be
   * given after the emission has been performed. The arguments are
   * the original \a dipole, the final dipole \a state and the emitted
   * \a parton. Also the scale (\a pt2) and \a type of the emission
   * must be supplied.
   */
  virtual bool
  finalVeto(tcEmiPtr dipole, Ariadne::tcDipoleStatePtr state, tcParPtr parton,
      Energy2 pt2, const EmissionType & type) const;

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
   * If positive, the phase space restriction on final-state emissions
   * is somewhat relaxed. If negative the phase space is strict and
   * only given by light-cone momentum non-ordering.
   */
  int isGenerous;

  /**
   * Veto only emissions from original partons.
   */
  bool useOnlyOriginal;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FSDipoleOrdering & operator=(const FSDipoleOrdering &);

};

}

#endif /* DIPSY_FSDipoleOrdering_H */
