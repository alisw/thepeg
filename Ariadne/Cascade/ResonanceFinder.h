// -*- C++ -*-
#ifndef Ariadne5_ResonanceFinder_H
#define Ariadne5_ResonanceFinder_H
//
// This is the declaration of the ResonanceFinder class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ResonanceFinder.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The ResonanceFinder class is used by the AriadneHandler to find
 * possible s-channel resonances in a given SubProcess. Normally this
 * is quite trivial, but in some cases the information about
 * resonances is not available in the SubProcess (eg. if the process
 * has been read in from an Les Houches event file). In this case the
 * user may inherit from this class to specify which of the final
 * state particles in the SubProcess was produced via a resonance.
 *
 * @see \ref ResonanceFinderInterfaces "The interfaces"
 * defined for ResonanceFinder.
 */
class ResonanceFinder: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ResonanceFinder();

  /**
   * The destructor.
   */
  virtual ~ResonanceFinder();
  //@}

public:

  /**
   * Return a list of intermediate s-channel resonances in the given
   * SubProcess. Sub-classes may insert resonances which they think
   * has been omitted in the event record. Resonances which have
   * decayed into further resonances must come before their children
   * in the list.
   */
  virtual tPVector resonances(SubProcess &) const;

  /**
   * Helper function to be used by sub-classes when inserting
   * resonances in the SubProcess. The \a resonance should be a newly
   * created particle and will be inserted in the \a sub process (and
   * the corresponding \a step). Mother/daughter relationships will be
   * modified. The \a children must already be included in the \a sub
   * process.
   */
  void insert(Step & step , SubProcess & sub,
	      PPtr resonance, const tPVector & children) const;

public:

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ResonanceFinder & operator=(const ResonanceFinder &);

};

}

#endif /* Ariadne5_ResonanceFinder_H */
