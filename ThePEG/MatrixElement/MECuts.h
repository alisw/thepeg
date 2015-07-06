// -*- C++ -*-
//
// MECuts.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_MECuts_H
#define ThePEG_MECuts_H
// This is the declaration of the MECuts class.

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Vectors/LorentzRotation.fh"
#include "ThePEG/Utilities/Triplet.h"
#include "ThePEG/PDT/StandardMatchers.h"

namespace ThePEG {

/**
 * The MECuts class is (currently not at all) used to make cuts on
 * generated phase space points for the hard interaction. A MECuts
 * object is selected for each matrix element. The EventHandler
 * has a default MECuts object, which may be overridden by the
 * selected SubProcessHandler object, which in turn may be overridden
 * by the selected MEBase object.
 *
 * The MECuts is used in two different ways. Individual handlers may
 * use the specific member functions which specify cuts on individual
 * variables. In addition the cut() member functions are always called
 * by the EventHandler to automatically check that all cuts are
 * passed. It is possible to derive new classes from the MECuts class,
 * in which case the virtual newcut() methods may be overridden and
 * will be called from the cut() methods.
 *
 *
 * @see EventHandler
 * @see SubProcessHandler
 * @see MEBase
 * @see Collision
 * @see SubProcess
 * 
 */
class MECuts: public Interfaced {

public:

  /**
   * Standard ctors and dtor
   */
  MECuts();

public:

  /**
   * This method is called by the EventHandler with the primary
   * SubProcess provided in its cm frame.
   */
  void cut(const SubProcess &) const;

public:

  /**
   * The minimum and maximum values of the invariant mass (squared) of
   * the hard sub-process.
   */
  Energy mHatMin() const { return theMHatMin; }
  /**
   * The minimum and maximum values of the invariant mass (squared) of
   * the hard sub-process.
   */
  Energy mHatMax() const { 
    theMHatMax > mHatMin()? theMHatMax: Constants::MaxEnergy; 
  }
  /**
   * The minimum and maximum values of the invariant mass (squared) of
   * the hard sub-process.
   */
  Energy sHatMin() const { return sqr(mHatMin()); }
  /**
   * The minimum and maximum values of the invariant mass (squared) of
   * the hard sub-process.
   */
  Energy sHatMax() const { return sqr(mHatMax()); }

  /**
   * The minimum and maximum values of the transverse momentum of the
   * outgoing particles in the hard sub-process.
   */
  Energy pTHatMin() const { return thePTHatMin; }
  /**
   * The minimum and maximum values of the transverse momentum of the
   * outgoing particles in the hard sub-process.
   */
  Energy pTHatMax() const { 
    return thePTHatMax > pTHatMin()? thePTHatMax: Constants::MaxEnergy; 
  }

  /**
   * Additional cut on the transverse momenta of the hard sub-process
   * for s-channel hard sub-processes for outgoing particles of mass
   * less than singularMassMax().
   */
  Energy pTHatSingularMin() const { return thePTHatSingularMin; }
  /**
   * Additional cut on the transverse momenta of the hard sub-process
   * for s-channel hard sub-processes for outgoing particles of mass
   * less than singularMassMax().
   */
  Energy singularMassMax() const { return theSingularMassMax; }

  /**
   * The minimum and maximum value of cosine of the scattering angle
   * in the restframe of a hard 2->2 scattering.
   */
  double cTHMin() const { return theCTHMin; }
  /**
   * The minimum and maximum value of cosine of the scattering angle
   * in the restframe of a hard 2->2 scattering.
   */
  double cTHMax() const { return theCTHMax; }

  /**
   * The minimum and maximum value of that of a hard 2->2 scattering.
   */
  Energy2 tHatMin() const { return theTHatMin; }
  /**
   * The minimum and maximum value of that of a hard 2->2 scattering.
   */
  Energy2 tHatMax() const { 
    return theTHatMax > tHatMin()? theTHatMax: Constants::MaxEnergy2; 
  }

  /**
   * The minimum and maximum value of uhat of a hard 2->2 scattering.
   */
  Energy2 uHatMin() const { return theUHatMin; }
  /**
   * The minimum and maximum value of uhat of a hard 2->2 scattering.
   */
  Energy2 uHatMax() const { 
    return theUHatMax > uHatMin()? theUHatMax: Constants::MaxEnergy2; 
  }

  /**
   * The minimum and maximum value of the scale in a hard scattering
   * as defined by the Handlers which performed the hard scattering.
   */
  Energy2 scaleMin() const { return theScaleMin; }
  /**
   * The minimum and maximum value of the scale in a hard scattering
   * as defined by the Handlers which performed the hard scattering.
   */
  Energy2 scaleMax() const {
    return theScaleMax > scaleMin()? theScaleMax: Constants::MaxEnergy2; 
  }

public:

  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

protected:

  /**
   * This method is called by the corresponding cut method with the
   * primary SubProcess provided in its cm frame. This bas-class
   * method does nothing.
   */
  virtual void newcut(const SubProcess &) const;

protected:

  /**
   * Standard Interfaced virtual functions.
   */
  virtual void doupdate();

  /**
   * Standard clone method.
   */
  virtual IBPtr clone() const;

private:

  /**
   * The minimum and maximum values of the invariant mass of
   * the hard sub-process.
   */
  Energy theMHatMin;
  /**
   * The minimum and maximum values of the invariant mass of
   * the hard sub-process.
   */
  Energy theMHatMax;

  /**
   * The minimum and maximum values of the transverse momentum of the
   * outgoing particles in the hard sub-process.
   */
  Energy thePTHatMin;
  /**
   * The minimum and maximum values of the transverse momentum of the
   * outgoing particles in the hard sub-process.
   */
  Energy thePTHatMax;

  /**
   * Additional cut on the transverse momenta of the hard sub-process
   * for s-channel hard sub-processes for outgoing particles of mass
   * less than theSingularMassMax.
   */
  Energy thePTHatSingularMin;
  /**
   * Additional cut on the transverse momenta of the hard sub-process
   * for s-channel hard sub-processes for outgoing particles of mass
   * less than theSingularMassMax.
   */
  Energy theSingularMassMax;

  /**
   * The minimum and maximum value of cosine of the scattering angle
   * in the restframe of a hard 2->2 scattering.
   */
  double theCTHMin;
  /**
   * The minimum and maximum value of cosine of the scattering angle
   * in the restframe of a hard 2->2 scattering.
   */
  double theCTHMax;

  /**
   * The minimum and maximum value of that of a hard 2->2 scattering.
   */
  Energy2 theTHatMin;
  /**
   * The minimum and maximum value of that of a hard 2->2 scattering.
   */
  Energy2 theTHatMax;

  /**
   * The minimum and maximum value of uhat of a hard 2->2 scattering.
   */
  Energy2 theUHatMin;
  /**
   * The minimum and maximum value of uhat of a hard 2->2 scattering.
   */
  Energy2 theUHatMax;

  /**
   * The minimum and maximum value of the scale in a hard scattering
   * as defined by the Handlers which performed the hard scattering.
   */
  Energy2 theScaleMin;
  /**
   * The minimum and maximum value of the scale in a hard scattering
   * as defined by the Handlers which performed the hard scattering.
   */
  Energy2 theScaleMax;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<MECuts> initMECuts;

  /**
   *  Private and non-existent assignment operator.
   */
  MECuts & operator=(const MECuts &);

};

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * MECuts.
 */
template <>
struct BaseClassTrait<MECuts,1>: public ClassTraitsType {
  /** Typedef of the base class of MECuts. */
  typedef Interfaced NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * MECuts class.
 */
template <>
struct ClassTraits<MECuts>:
  /** Return the class name. */
    public ClassTraitsBase<MECuts> {
  static string className() { return "ThePEG::MECuts"; }
};

/** @endcond */

}

#endif /* ThePEG_MECuts_H */
