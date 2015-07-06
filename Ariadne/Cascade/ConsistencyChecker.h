// -*- C++ -*-
#ifndef Ariadne5_ConsistencyChecker_H
#define Ariadne5_ConsistencyChecker_H
//
// This is the declaration of the ConsistencyChecker class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "DipoleState.fh"
#include "Emission.fh"
#include "Parton.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The ConsistencyChecker class is used by the AriadneHandler to check
 * a DipoleState for consistency. The check is done on the initial
 * DipoleState, and then after each Emission, to veto any emission
 * which gives an inconsistent state. If the initial DipoleState is
 * inconsistent, the whole event may be vetoed, otherwise further
 * checking is suspended. This base class optionally checks that all
 * gluons have an invariant transverse momentum above the cutoff and
 * that all stings are large enough to be able to produce a minimum
 * set of hadrons. Sub-classes may introduce further checks.
 *
 * @see \ref ConsistencyCheckerInterfaces "The interfaces"
 * defined for ConsistencyChecker.
 */
class ConsistencyChecker: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ConsistencyChecker(): vetoInitial(true), tryFixing(false), checkGluonPT(false),
			checkStrings(false) {}

  /**
   * The destructor.
   */
  virtual ~ConsistencyChecker();
  //@}

protected:

  /**
   * The main virtual function. Returns true if the given DipoleState
   * is consistent. The initial DipoleState is called with the
   * Emission pointer set to null.
   */
  virtual bool isConsistent(const DipoleState &, tcEmPtr) const;

  /**
   * If the initial state is inconsistent but not vetoed, this
   * function may be implemented to possibly fix the state to become
   * consistent. Return false if it could not be fixed.
   */
  virtual bool fixme(DipoleState &) const {
    return false;
  }

public:

  /**
   * The function which is called by the AriadneHandler. Forwards the
   * call to isConsistent() and, in case of the initial DipoleState,
   * optionally vetoes the whole event.
   */
  bool check(DipoleState & state, tcEmPtr emission) const {
    if ( isConsistent(state, emission) ) return true;
    if ( !emission && tryFixing && fixme(state) ) return true;
    if ( !emission && vetoInitial ) throw Veto();
    return false;
  }

protected:

  /** 
   * Returns true if all gluon have an invariant transverse momentum
   * above the cutoff.
   */
  bool gluonPT(const DipoleState &) const;

  /** 
   * Returns false if a gluon \a p2 has an invariant transverse
   * momentum below \a pt2cut w.r.t. \a p1 and \a p3. If none of the
   * partons have been touched or any of the partons are failsafe, no
   * check is done and true is returned.
   */
  bool gluonPT(tcParPtr p1, tcParPtr p2, tcParPtr p3, Energy2 pt2cut) const;

  /**
   * Returns true if all remnants are consistent.
   */
  bool remnantPhaseSpace(const DipoleState & state) const;

  /**
   * Returns true if all strings have large enough mass.
   */
  bool strings(const DipoleState &) const;

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
   * If true, an event is vetoed if the initial DipoleState is
   * inconsistent.
   */
  bool vetoInitial;

  /**
   * If true, a broken state may be fixed by subclasses.
   */
  bool tryFixing;

  /**
   * If true, all gluons are required to have an invariant transverse
   * momentum above the cutoff.
   */
  bool checkGluonPT;

  /**
   * If true, all strings are required to have an invariant mass big
   * enough to produce a minimum set of hadrons.
   */
  bool checkStrings;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ConsistencyChecker & operator=(const ConsistencyChecker &);

};

}

#endif /* Ariadne5_ConsistencyChecker_H */
