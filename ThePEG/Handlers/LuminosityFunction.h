// -*- C++ -*-
//
// LuminosityFunction.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_LuminosityFunction_H
#define ThePEG_LuminosityFunction_H
// This is the declaration of the LuminosityFunction class.

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/LastXCombInfo.h"
#include "ThePEG/Vectors/LorentzRotation.fh"
#include "ThePEG/Utilities/Interval.h"

namespace ThePEG {

/**
 * The LuminosityFunction describes the momentum distribution of the
 * incoming beams in an experiment. This is used by a EventHandler to
 * generate collisions in their CM system. The LuminosityFunction will
 * be asked to produce a LorentzRotation giving the transformation to
 * the laboratory system.
 *
 * The LuminosityFunction inherits from the LastXCombInfo class to
 * give easy access to the information of the generated primary
 * sub-process in the selected XComb.
 *
 * This base class implements simple fixed momentum beams with
 * energies given by the BeamEMaxA and BeamEMaxB interfaces.
 *
 * @see \ref LuminosityFunctionInterfaces "The interfaces"
 * defined for LuminosityFunction.
 * @see XComb
 * 
 */
class LuminosityFunction: public HandlerBase, public LastXCombInfo<> {

  /** EventHandler is a friend. */
  friend class EventHandler;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor. Optionally the maximum energy of beam \a a
   * and \a b can be given.
   */
  LuminosityFunction(Energy a = 45.6*GeV, Energy b = 45.6*GeV);
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

  /**
   * How many random numbers are needed to generate a phase space
   * point? Default is zero in which means the energy of the incoming
   * particles is fixed. The only other reasonable values are 1 and 2.
   */
  virtual int nDim(const cPDPair &) const;

  /**
   * The value of the luminosity function for the given particle types
   * for the given energy fractions l1 and l2 (\f$l=\log(1/x)\f$). The
   * default version returns 1 if l1 and l2 are zero otherwize zero.
   */
  virtual double value(const cPDPair &, double l1, double l2) const;

  /**
   * Generate energy fractions l1 and l2 (\f$l=\log(1/x)\f$) given
   * 'nDim()' random numbers in the range ]0,1[ given by the
   * iterators. The jacobian argument must be multiplied by the
   * jacobian of the variable transformation to l1 and l2. The default
   * version is just a delta function with a jacobian of 1.
   */
  virtual pair<double,double>
  generateLL(const double * r, double & jacobian) const;
  //@}

public:

  /** @name Simple access functions */
  //@{
  /**
   * The maximum energy of the beam entering along the positive z-axis.
   */
  Energy beamEMaxA() const { return theBeamEMaxA; }

  /**
   * The maximum energy of the beam entering along the negative z-axis.
   */
  Energy beamEMaxB() const { return theBeamEMaxB; }
  //@}

protected:

  /**
   * The maximum energy of the beam entering along the positive z-axis.
   */
  void beamEMaxA(Energy x) { theBeamEMaxA = x; }

  /**
   * The maximum energy of the beam entering along the negative z-axis.
   */
  void beamEMaxB(Energy x) { theBeamEMaxB = x; }

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
   * Standard Init function used to initialize the interface.
   */
  static void Init();

  /**
   * Set information about the selected XComb.
   */
  void select(tXCombPtr);

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
   * The maximum energy of the beam entering along the positive z-axis.
   */
  Energy theBeamEMaxA;

  /**
   * The maximum energy of the beam entering along the negative z-axis.
   */
  Energy theBeamEMaxB;

private:

  /**
   * Describe an abstract class with persistent data.
   */
  static ClassDescription<LuminosityFunction> initLuminosityFunction;

  /**
   *  Private and non-existent assignment operator.
   */
  LuminosityFunction & operator=(const LuminosityFunction &);

};

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of LuminosityFunction.
 */
template <>
struct BaseClassTrait<LuminosityFunction,1>: public ClassTraitsType {
  /** Typedef of the base class of LuminosityFunction. */
  typedef HandlerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * LuminosityFunction class.
 */
template <>
struct ClassTraits<LuminosityFunction>:
    public ClassTraitsBase<LuminosityFunction> {
  /** Return the class name. */
  static string className() { return "ThePEG::LuminosityFunction"; }
};

/** @endcond */

}

#endif /* ThePEG_LuminosityFunction_H */
