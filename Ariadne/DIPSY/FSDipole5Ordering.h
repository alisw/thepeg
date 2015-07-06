// -*- C++ -*-
#ifndef DIPSY_FSDipole5Ordering_H
#define DIPSY_FSDipole5Ordering_H
//
// This is the declaration of the FSDipole5Ordering class.
//

#include "Ariadne/Cascade/ReweightBase.h"

namespace DIPSY {

using namespace Ariadne5;

/**
 * Here is the documentation of the FSDipole5Ordering class.
 *
 * @see \ref FSDipole5OrderingInterfaces "The interfaces"
 * defined for FSDipole5Ordering.
 */
class FSDipole5Ordering: public Ariadne5::ReweightBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FSDipole5Ordering();

  /**
   * The destructor.
   */
  virtual ~FSDipole5Ordering();
  //@}

public:

  /**
   * In addition to the reweight function a final hit/miss veto may be
   * given after the \a emission has been performed. Will only be
   * called if hasFinalVeto() returns true.
   */
  virtual bool finalVeto(const Emission & emission) const;

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

  /**
   * Fudge factor for +/- ordering. Values above one increases the
   * phase space.
   */
  double f;

  /**
   * Limit on invariant transverse momentum. Any emission below this
   * is not checked.
   */
  Energy ptmin;

  /**
   * Option for suppression of hard radiation if partons in the dipole
   * are far separated.
   */
  int hardSuppression;

  /**
   * Fudge factor to scale the fraction of light-cone momentum
   * available for radiation.
   */
  double fudge;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FSDipole5Ordering & operator=(const FSDipole5Ordering &);

};

}

#endif /* DIPSY_FSDipoleOrdering_H */
