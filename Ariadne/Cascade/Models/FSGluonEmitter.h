// -*- C++ -*-
#ifndef ARIADNE5_FSGluonEmitter_H
#define ARIADNE5_FSGluonEmitter_H
//
// This is the declaration of the FSGluonEmitter class.
//

#include "Ariadne/Cascade/EmitterBase.h"
#include "FSGluonEmission.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The FSGluonEmitter class implements the standard classical gluon
 * emission from a final-state colour dipole.
 *
 * @see \ref FSGluonEmitterInterfaces "The interfaces"
 * defined for FSGluonEmitter.
 */
class FSGluonEmitter: public EmitterBase {

public:

  /**
   * Convenient typedef.
   */
  ThePEG_DECLARE_POINTERS(Ariadne5::FSGluonEmission,FSGEmPtr);

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FSGluonEmitter();

  /**
   * The destructor.
   */
  virtual ~FSGluonEmitter();
  //@}

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return true if and only if this model can handle the given
   * Emitter.
   */
  virtual bool canHandle(const DipoleBase &) const;

  /**
   * If the given EmissionModel overlaps with this model for the given
   * Emitter, return true if this model should take precedence. Must
   * only be called for Emitters for which canHandle() is true.
   */
  virtual bool overrides(const EmitterBase &, DipoleBase &) const;

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
   * Generic function to generate a phase space point for a gluon
   * emission, \a e, given a maximum and minimum evolution scale, the
   * total availabel energy, \a W, a prefactor, \a C, the maximum
   * rapidity interval, \a yint, and the power of the exponents in the
   * dipole splitting function.
   *
   * @return a negative value if no emission was generated in the
   * given interval; zero if a point has been generated but should
   * subsequently be thrown aaway; and a number between 0 and one if a
   * point has been generated which should be kept with the
   * corresponding probability.
   */
  static double gluonEngine(FSGluonEmission & e, Energy rhomin, Energy rhomax,
			    Energy W, double C, double yint, int n1, int n3);

  /**
   * Generic function for inserting a gluon in a QCDDipole, \a d, this
   * creating a new dipole. The parton \a p will be assumed to be the
   * emitting one and the gluon will inherit its corresponding (anti-)
   * colour.
   *
   * @return the created gluon and the created QCDDipole.
   */
  static pair<ParPtr, QCDPtr> insertGluon(QCDDipole & d, tParPtr p);

  /**
   * Remove a gluon previously added by insertGluon() 
   */
  static void removeGluon(QCDDipole & d, tParPtr g, tParPtr p);

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
  FSGluonEmitter & operator=(const FSGluonEmitter &);

};

}

#endif /* ARIADNE5_FSGluonEmitter_H */
