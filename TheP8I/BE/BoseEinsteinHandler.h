// -*- C++ -*-
#ifndef THEP8I_BoseEinsteinHandler_H
#define THEP8I_BoseEinsteinHandler_H
//
// This is the declaration of the BoseEinsteinHandler class.
//

#include "ThePEG/Handlers/DecayHandler.h"
#include "TheP8I/Config/Pythia8Interface.h"
#include "BoseEinstein.h"

namespace TheP8I {

using namespace ThePEG;

/**
 * The BoseEinsteinHandler inherits from the basic DecayHandler class
 * and divides the decay process into three steps. First the
 * short-lived particles are decayed (as defined by
 * <interface>BoseEinstein_widthSep</interface>), then all particle
 * momenta are shifted to simulate Bose-Einstein correlation effects,
 * whereafter the decay chains are continued as in a normal
 * DecayHandler.
 *
 * Note that if this class is used as a PreDecayHandler, it may be
 * necessary to switch off the <interface>ContinueDecays</interface>
 * flag.
 *
 * @see \ref BoseEinsteinHandlerInterfaces "The interfaces"
 * defined for BoseEinsteinHandler.
 */
class BoseEinsteinHandler: public DecayHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  BoseEinsteinHandler();

  /**
   * The destructor.
   */
  virtual ~BoseEinsteinHandler();
  //@}

public:

  /**
    * The main function called by the EventHandler class to
    * perform the Decay step.
    * @param eh the EventHandler in charge of the Event generation.
    * @param tagged if not empty these are the only particles which should
    * be considered by the StepHandler.
    * @param hint a Hint object with possible information from previously
    * performed steps.
    * @throws Veto if the StepHandler requires the current step to be
    * discarded.
    * @throws Stop if the generation of the current Event should be stopped
    * after this call.
    * @throws Exception if something goes wrong.
    */
  virtual void handle(EventHandler & eh, const tPVector & tagged,
		      const Hint & hint);

  /**
   * Perform the decay of one unstable particle if its width is large enough.
   * @param parent the particle to be decayed.
   * @param s the Step where decay products are inserted.
   */
  void performDecay(tPPtr parent, Step & s) const;

  /**
   * Add a particle to the Pythia8::Event. If the particle has
   * decayed, instead add its children (recursively).
   */
  void addParticle(tPPtr p);

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



protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * The interface to the Pythia 8 object,
   */
  Pythia8Interface pythia;

  /**
   * The BoseEinstein object from Pythia8.
   */
  Pythia8::BoseEinstein theBE;

  /**
   * The particles sent to Pythia8 for BE-shifting
   */
  tPVector particles;

  /**
   * Continue the decay chains after the Bose-Einstein shifting.
   */
  bool doContinueDecays;

#include "BoseEinsteinHandler-var.h"

public:

  /**
   * Exception class declaration.
   */
  struct BoseEinsteinError: public Exception {};

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<BoseEinsteinHandler> initBoseEinsteinHandler;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BoseEinsteinHandler & operator=(const BoseEinsteinHandler &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of BoseEinsteinHandler. */
template <>
struct BaseClassTrait<TheP8I::BoseEinsteinHandler,1> {
  /** Typedef of the first base class of BoseEinsteinHandler. */
  typedef DecayHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BoseEinsteinHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<TheP8I::BoseEinsteinHandler>
  : public ClassTraitsBase<TheP8I::BoseEinsteinHandler> {
  /** Return a platform-independent class name */
  static string className() { return "TheP8I::BoseEinsteinHandler"; }
  /**
   * The name of a file containing the dynamic library where the class
   * BoseEinsteinHandler is implemented. It may also include several, space-separated,
   * libraries if the class BoseEinsteinHandler depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libTheP8I.so"; }
};

/** @endcond */

}

#endif /* THEP8I_BoseEinsteinHandler_H */
