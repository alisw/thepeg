// -*- C++ -*-
//
// FixedCMSLuminosity.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_FixedCMSLuminosity_H
#define ThePEG_FixedCMSLuminosity_H
// This is the declaration of the FixedCMSLuminosity class.

#include "LuminosityFunction.h"

namespace ThePEG {

/**
 * The FixedCMSLuminosity class describes an experiment with incoming
 * particles colliding with precicely defined and opposite momenta. It
 * is derived from the LuminosityFunction base class.
 *
 * \deprecated As the LuminosityFunction base class has increased
 * functionality (exceeding the functionality of this class) the use
 * of FixedCMSLuminosity is deprecated, and the class will be removed
 * in a future release. Note also that by setting the individual beam
 * energies in the base class, the behavior of this object may be
 * inconsistent, in that the collision will not, as specified, be in
 * the center-of-mass system.
 *
 * @see \ref FixedCMSLuminosityInterfaces "The interfaces"
 * defined for FixedCMSLuminosity.
 */
class FixedCMSLuminosity: public LuminosityFunction {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Destructor.
   */
  virtual ~FixedCMSLuminosity();
  //@}

public:

  /**
   * The total energy in the cms of the incoming particles.
   */
  Energy energy() const { return maximumCMEnergy(); }

public:

  /**
   * Standard Init function used to initialize the interface.
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

private:

  /**
   * Utility function used by the interface.
   */
  void setEnergy(Energy);

  /**
   * Utility function used by the interface.
   */
  Energy getEnergy() const;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static NoPIOClassDescription<FixedCMSLuminosity> initFixedCMSLuminosity;

  /**
   *  Private and non-existent assignment operator.
   */
  FixedCMSLuminosity & operator=(const FixedCMSLuminosity &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of FixedCMSLuminosity.
 */
template <>
struct BaseClassTrait<FixedCMSLuminosity,1>: public ClassTraitsType {
  /** Typedef of the base class of FixedCMSLuminosity. */
  typedef LuminosityFunction NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * FixedCMSLuminosity class and the shared object where it is defined.
 */
template <>
struct ClassTraits<FixedCMSLuminosity>:
    public ClassTraitsBase<FixedCMSLuminosity> {
  /** Return the class name. */
  static string className() { return "ThePEG::FixedCMSLuminosity"; }
  /** Return the name of the shared library be loaded to get access to
   *  the FixedCMSLuminosity class and every other class it uses
   *  (except the base class). */
  static string library() { return "FixedCMSLuminosity.so"; }
};

/** @endcond */

}

#endif /* ThePEG_FixedCMSLuminosity_H */
