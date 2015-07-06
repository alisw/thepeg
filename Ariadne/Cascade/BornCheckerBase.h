// -*- C++ -*-
#ifndef Ariadne5_BornCheckerBase_H
#define Ariadne5_BornCheckerBase_H
//
// This is the declaration of the BornCheckerBase class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "BornCheckerBase.fh"
#include "DipoleState.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The BornCheckerBase is the base class for any class implementing
 * the checking of a DipoleState reclustered in the CKKW-L algorithm,
 * to make sure it corresponds to a reasonable Born-level state (the
 * lowest jet-multiplicity state in a CKKW-L merging).
 *
 * @see \ref BornCheckerBaseInterfaces "The interfaces"
 * defined for BornCheckerBase.
 */
class BornCheckerBase: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  BornCheckerBase();

  /**
   * The destructor.
   */
  virtual ~BornCheckerBase();
  //@}

public:

  /** @name Virtual functions to be overidden in sub-classes. */
  //@{
  /**
   * If the given state corresponds to a valid Born-level, minimum
   * jet-multiplicity in a CKKW-L merging, return the factorization
   * scale associated with that state. If the given state does not
   * correspond to a Born-level state, return the the factorization
   * scale, but with a minus sign. If the given state is completely
   * unknown return zero. The factorization scale should be the one
   * normally used when Ariadne is run non-merged.
   */
  virtual Energy check(const DipoleState &) const = 0;

  /**
   * If the given state corresponds to a valid state from a matrix
   * element generator in a CKKW-L merging return true.
   */
  virtual bool checkTree(const DipoleState &) const = 0;

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
  BornCheckerBase & operator=(const BornCheckerBase &);

};

}

#endif /* Ariadne5_BornCheckerBase_H */
