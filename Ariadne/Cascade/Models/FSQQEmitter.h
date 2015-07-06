// -*- C++ -*-
#ifndef ARIADNE5_FSQQEmitter_H
#define ARIADNE5_FSQQEmitter_H
//
// This is the declaration of the FSQQEmitter class.
//

#include "Ariadne/Cascade/EmitterBase.h"
#include "FSQQEmission.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The FSQQEmitter class implements the final-state splitting of gluons
 * into a q-qbar pair.
 *
 * @see \ref FSQQEmitterInterfaces "The interfaces"
 * defined for FSQQEmitter.
 */
class FSQQEmitter: public EmitterBase {

public:

  /**
   * Convenient typedef.
   */
  ThePEG_DECLARE_POINTERS(Ariadne5::FSQQEmission,FSQQEmPtr);

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FSQQEmitter();

  /**
   * The destructor.
   */
  virtual ~FSQQEmitter();
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
   * Check if objects related to the given \a dipole have been touched
   * in a way such that emissions must be regenerated.
   */
  virtual bool touched(const DipoleBase & dipole) const;

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
   * splitting, \a e, given a maximum and minimum evolution scale, the
   * total availabel energy, \a W, a prefactor, \a C, the maximum
   * rapidity interval, \a yint, and the scaled masses of the partons,
   * \a sy1 and \a syq.
   *
   * @return a negative value if no emission was generated in the
   * given interval; zero if a point has been generated but should
   * subsequently be thrown aaway; and a number between 0 and one if a
   * point has been generated which should be kept with the
   * corresponding probability.
   */
  static double
  qqbarEngine(FSQQEmission & e, Energy rhomin, Energy rhomax,
	      Energy W, double C, double yint, double sy1, double syq);

  /**
   * Split the given \a gluon in a \a dipole into a q-qbar pair with
   * flavour ifl.
   *
   * @return a pair of partons where the first is the quark and the
   * second is the anti-quark.
   */
  static pair<tParPtr,tParPtr>
  splitGluon(QCDDipole & dipole, tParPtr gluon, int ifl);

  /**
   * Fuse a \a q and \a qbar into a gluon as if previously split by
   * splitGluon. If actually a previous split, \a gluon must be
   * non-zero. \a d is the dipole connecting to either the \a q or the
   * \a qbar, while \a od is the other dipole.
   *
   * @return the gluon.
   */
  static tParPtr
  fuseQQBar(QCDDipole & d, tParPtr q, tParPtr qbar, tQCDPtr od, tParPtr gluon);

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
  FSQQEmitter & operator=(const FSQQEmitter &);

};

}

#endif /* ARIADNE5_FSQQEmitter_H */
