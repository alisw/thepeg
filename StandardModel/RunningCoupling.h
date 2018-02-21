// -*- C++ -*-
//
// RunningCoupling.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad, (C) 2009 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_RunningCoupling_H
#define ThePEG_RunningCoupling_H
// This is the declaration of the RunningCoupling class.

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "StandardModelBase.fh"

namespace ThePEG {

/**
 * RunningCoupling an abstract base class unifying the treatment
 * of running couplings in ThePEG.
 *
 * @see \ref RunningCouplingInterfaces "The interfaces"
 * defined for RunningCoupling.
 * @see StandardModelBase
 */
class RunningCoupling: public Interfaced {

public:

  /**
   * The default constructor.
   */
  RunningCoupling () : theScaleFactor(1.) {}

  /**@name Methods to be implemented by a derived class */
  //@{

  /**
   * Return the value of the coupling at a given \a scale using the
   * given standard model object, \a sm.
   */
  virtual double value (Energy2 scale, const StandardModelBase & sm) const = 0;

  /**
   * Return the number of loops contributing to
   * the running this coupling. The default returns
   * zero to ensure backward compatibility.
   */
  virtual unsigned int nloops () const { return 0; }

  //@}

  /**
   * Return the value of the coupling at a given \a scale using the
   * StandardModelBase object used by the EventGenerator.
   */
  double value(Energy2 scale) const {
    return value(scale,*(generator()->standardModel()));
  }

  /**
   * Return an overestimate to the running coupling at the
   * given scale. This is defined to aid veto algorithms
   * and by default returns the coupling itself, using the EventGenerators's
   * StandardModelBase object.
   */
  virtual double overestimateValue (Energy2 scale) const {
    return value(scale);
  }

  /**
   * Return the ratio of the exact to the overestimated value
   * of the running coupling. The default implementation returns
   * one in accordance with the default implementation of
   * overestimateValue
   */
  virtual double ratioToOverestimate (Energy2) const {
    return 1.;
  }

  /**
   * Return the scale factor, which may be used to globally
   * rescale the argument of the running coupling.
   */
  double scaleFactor () const { return theScaleFactor; }

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

private:

  /**
   * Describe an abstract class without persistent data.
   */
  static AbstractClassDescription<RunningCoupling> initRunningCoupling;

  /**
   *  Private and non-existent assignment operator.
   */
  RunningCoupling & operator=(const RunningCoupling &);

  /**
   * The scale factor used to rescale the argument of
   * the running coupling.
   */
  double theScaleFactor;

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of RunningCoupling. */
template <>
struct BaseClassTrait<RunningCoupling,1>: public ClassTraitsType {
  /** Typedef of the first base class of RunningCoupling. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  RunningCoupling class. */
template <>
struct ClassTraits<RunningCoupling>: public ClassTraitsBase<RunningCoupling> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::RunningCoupling"; }
};

/** @endcond */

}

#endif /* ThePEG_RunningCoupling_H */
