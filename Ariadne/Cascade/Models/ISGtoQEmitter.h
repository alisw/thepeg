// -*- C++ -*-
#ifndef ARIADNE5_ISGtoQEmitter_H
#define ARIADNE5_ISGtoQEmitter_H
//
// This is the declaration of the ISGtoQEmitter class.
//

#include "ISQEmitter.h"
#include "ISGtoQEmission.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The ISGtoQEmitter class implements the initial-state splitting of gluons
 * into a q-qbar pair.
 *
 * @see \ref ISGtoQEmitterInterfaces "The interfaces"
 * defined for ISGtoQEmitter.
 */
class ISGtoQEmitter: public ISQEmitter {

public:

  /**
   * Convenient typedef.
   */
  ThePEG_DECLARE_POINTERS(Ariadne5::ISGtoQEmission,ISGtoQEmPtr);

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ISGtoQEmitter();

  /**
   * The destructor.
   */
  virtual ~ISGtoQEmitter();
  //@}

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return true if and only if this model can handle the given
   * Emitter.
   */
  virtual bool canHandle(const DipoleBase &) const;

  /**
   * Generate the a phase space point for an emission corresponding to
   * this model. Must only be called for Emitters for which
   * canHandle() is true.
   */
  virtual EmPtr generate(const DipoleBase &, Energy rhomin, Energy rhomax) const;

  /**
   * Perform an emission previously generated for this Emitter. Must
   * only be called for Emitters for which canHandle() is true.
   * @return true if the emission was successful
   */
  virtual bool perform(const Emission &) const;

  /**
   * Reverse a previously performed emission. Sub-classes which has
   * signalled that they can revert an emission but fails to do so,
   * must throw a Exception::runerror.
   */
  virtual void revert(const Emission & emission) const;

  //@}

public:

  /**
   * Helper function for the main generate function implementing
   * genetation from a given remnant parton.
   */
  virtual EmPtr generate(tRemParPtr rem, const QCDDipole & d,
			       Energy rhomin, Energy rhomax) const;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ISGtoQEmitter & operator=(const ISGtoQEmitter &);

};

}

#endif /* ARIADNE5_ISGtoQEmitter_H */
