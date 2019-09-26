// -*- C++ -*-
#ifndef ThePEG_FixedTargetLuminosity_H
#define ThePEG_FixedTargetLuminosity_H
//
// This is the declaration of the FixedTargetLuminosity class.
//

#include "LuminosityFunction.h"

namespace ThePEG {

/**
 * Here is the documentation of the FixedTargetLuminosity class.
 *
 * @see \ref FixedTargetLuminosityInterfaces "The interfaces"
 * defined for FixedTargetLuminosity.
 */
class FixedTargetLuminosity: public LuminosityFunction {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FixedTargetLuminosity();

  /**
   * The destructor.
   */
  virtual ~FixedTargetLuminosity();
  //@}
  
  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return true if this luminosity function can actually handle a
   * given pair of incoming particles.
   */
  virtual bool canHandle(const cPDPair &) const;

  /**
   * Return the maximum possible center of mass energy for an event.
   */
  virtual Energy maximumCMEnergy() const;

  /**
   * Return the rotation needed to transform from the collision cm
   * system to the labotatory system. This default version returns the
   * unit transformation.
   */
  virtual LorentzRotation getBoost() const;

  /**
   * Return the rapidity of the colliding particles (at the maximum
   * energy) in the laboratory system. This default version assumes
   * the CM system is the same as the lab system and returns zero.
   */
  virtual double Y() const;
  //@}
  
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
protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}
  
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FixedTargetLuminosity & operator=(const FixedTargetLuminosity &) = delete;

private:

  /**
   *    The beam particle
   */
  PDPtr beam_;

  /**
   *    The target particle
   */
  PDPtr target_;

  /**
   *  CMS energy
   */
  Energy ecms_;

  /**
   *  Boost
   */
  double beta_;
  
  
};

}

#endif /* ThePEG_FixedTargetLuminosity_H */
