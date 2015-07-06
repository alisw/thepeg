// -*- C++ -*-
#ifndef PYTHIA7_ShowerHandler_H
#define PYTHIA7_ShowerHandler_H
// This is the declaration of the ShowerHandler class.

#include "Pythia7/Config/Pythia7.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/Utilities/ObjectIndexer.h"
// #include "ShowerHandler.fh"
// #include "ShowerHandler.xh"

#include "Pythia7/Shower/SpaceShowerHandler.h"
#include "Pythia7/Shower/TimeShowerHandler.h"
#include "Pythia7/Shower/Basics.h"
#include "Pythia7/Shower/Shower.h"
#include "ThePEG/PDF/PartonBinInstance.h"

namespace Pythia7 {

/**
 * The <code>ShowerHandler</code> class administers parton showers for
 * external partons in a hard SubProcess or in a
 * decay of a heavy resonance. The main function will identify partons
 * in the current step which should shower and create corresponding
 * TimeParticle and SpaceParticle
 * objects which, using the assigned TimeShower and
 * SpaceShower objects, will administer the
 * showering.
 *
 * @see \ref ShowerHandlerInterfaces "The interfaces"
 * defined for ShowerHandler.
 * @see SubProcess
 * @see TimeParticle
 * @see SpaceParticle
 * @see TimeShower
 * @see SpaceShower
 * 
 */
class ShowerHandler: public CascadeHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ShowerHandler();

  /**
   * The copy constructor.
   */
  inline ShowerHandler(const ShowerHandler &);

  /**
   * The destructor.
   */
  virtual ~ShowerHandler();
  //@}

public:

  /** @name Virtual functions required by the CascadeHandler class. */
  //@{
  /**
    * The main function called by the EventHandler class to
    * perform a step.  Stores important information and calls
    * cascade().
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
   * Perform the actual cascade.
   */
  virtual void cascade();
  //@}

  /**
   * Perform a final-state cascade on the given partons.
   */
  virtual void cascade(const tPVector & final);

  /**
   * Perform an initial-state cascade on the incoming partons of the
   * given sub-process.
   */
  virtual void cascade(tSubProPtr sub);

  /**
   * Return true if the given particle can be showered. If \a inc is
   * true, test if the particle can initiate an initial state shower.
   */
  virtual bool canShower(tcPPtr p, bool inc = false) const;

  /**
   * If true all branchings in the shower will be added to the new
   * step. Otherwise only the final partons will be added.
   */
  inline bool history() const;

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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /**
   * Add a particle to the internal event record with the given
   * status. Note that if the parent particles has not already been
   * added, they will not be listed as parents in the internal event
   * record.
   */
  long addParticle(tPPtr p, long status);
  /**
   * Add a particle to the internal event record with the given status
   * and mother indices. Note that if the parent particles has not
   * already been added, they will not be listed as parents in the
   * internal event record.
   */
  long addParticle(tPPtr p, long status, long mother1, long mother2);

  /**
   * Add a time-like particle to the internal event record with the
   * given mother indices. If the particle has decayed add its
   * children recursively.
   */
  void addFinalParticles(tPPtr p, long mother1, long mother2);

  /**
   * Set the static parameters of the underlying model object and
   * return an instance.
   */
  Shower::Shower * getModel();

  /**
   * Return true if the given particle is a resonance (defined by the
   * fact that all children only has the resonance as a parent).
   */
  bool isResonance(tcPPtr r) const;

protected:

  /**
   * Recursively find the first and last indices in the internal event
   * record of the incoming shower initiator.
   */
  void findChain(pair<long,long> & chain, long init) const;

  /**
   * Find the parent ThePEG::Particle pair of a given entry in the
   * internal event record.
   */
  tPPair findParent(long i) const;

  /**
   * Create a ThePEG::Particle corresponding to a given entry in the
   * internal event record. If \a inc is true the particle is an
   * initiator for an intial state shower.
   */
  tPPtr createParticle(long i, bool inc = false);

  /**
   * Copy information of the given entry in the internal event record
   * to the given ThePEG::Particle.
   */
  tPPtr copyParticle(long i, tPPtr p);

  /**
   * Return the ThePEG::Particle corresponding to a given entry in the
   * internal event record.
   */
  tPPtr getParticle(long i) const;

  /**
   * Return the momentum of the given entry in the internal event
   * record.
   */
  Lorentz5Momentum momentum(long i) const;

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

public:

  /**
   * Return a pointer to the object which implements the actual
   * time-like shower model.
   */
  inline tTimeShowerPtr timeShower() const;
  /**
   * Return a pointer to the object which implements the actual
   * space-like shower model.
   */
  inline tSpaceShowerPtr spaceShower() const;

protected:


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * The object implementing the administration of the time- and
   * space-like showers.
   */
  Shower::Shower * theShowerModel;

  /**
   * The object which implements the actual time-like shower model.
   */
  TimeShowerPtr theTimeShower;

  /**
   * The object which implements the actual space-like shower model.
   */
  SpaceShowerPtr theSpaceShower;

  /**
   * If true add final-state radiation off initial-state cascade.
   */
  bool addFSROnISR;

  /**
   * If true all branchings in the shower will be added to the new
   * step. Otherwise only the final partons will be added.
   */
  bool theFullHistory;

  /**
   * The maximum number of attempts to cascade a sub process before
   * giving up and throwing an eventerror.
   */
  long maxTries;

protected:

  /**
   * The event record used internally.
   */
  Shower::Event event;

  /**
   * Association between ColourLines and colour indices.
   */
  ObjectIndexer<long,ColourLine> colourIndex;

  /**
   * Association between Particles and indices.
   */
  ObjectIndexer<long,Particle> particleIndex;

private:

  /**
   * Exception class used if too many attempts is made to shower a
   * sub-process.
   */
  struct InfiniteLoopException: public Exception {};
private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ShowerHandler> initShowerHandler;

  /**
   *  Private and non-existent assignment operator.
   */
  ShowerHandler & operator=(const ShowerHandler &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * Pythia7::ShowerHandler.
 */
template <>
struct BaseClassTrait<Pythia7::ShowerHandler,1>: public ClassTraitsType {
  /** Typedef of the base class of Pythia7::ShowerHandler. */
  typedef CascadeHandler NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Pythia7::ShowerHandler class and the shared object where it is
 * defined.
 */
template <>
struct ClassTraits<Pythia7::ShowerHandler>
  : public ClassTraitsBase<Pythia7::ShowerHandler> {
  /** Return the class name. */
  static string className() { return "Pythia7::ShowerHandler"; }
  /** Return the name of the shared library be loaded to get access to
   *  the Pythia7::ShowerHandler class and every other class it uses
   *  (except the base class). */
  static string library() { return "libP7Shower.so"; }

};

/** @endcond */

}

#include "ShowerHandler.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "ShowerHandler.tcc"
#endif

#endif /* PYTHIA7_ShowerHandler_H */
