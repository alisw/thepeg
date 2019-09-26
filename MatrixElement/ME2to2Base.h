// -*- C++ -*-
//
// ME2to2Base.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ME2to2Base_H
#define ThePEG_ME2to2Base_H
// This is the declaration of the ME2to2Base class.

#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Interface/Switch.fh"

namespace ThePEG {

/**
 * ME2to2Base can be used as a base class for any matrix element class
 * implementing 2\f$\rightarrow\f$ 2 processes. It extends the MEBase
 * base class by implementing the virtual scale() method to return the
 * assumed scale of a given process according to several options. It
 * also caches some useful characteristics of the chosen phase space
 * point, such as tHat(), uHat() and the masses of the external
 * partons.
 *
 * @see \ref ME2to2BaseInterfaces "The interfaces"
 * defined for ME2to2Base.
 * @see MEBase
 */
class ME2to2Base: public MEBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  ME2to2Base()
    : theScaleChoice(0), theLastTHat(ZERO), theLastUHat(ZERO),
      theLastPhi(0.0) {}

  /**
   * Destructor.
   */
  virtual ~ME2to2Base();
  //@}

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * The number of internal degreed of freedom used in the matrix
   * element.
   */
  virtual int nDim() const { return 1; }

  /**
   * Generate internal degrees of freedom given 'nDim()' uniform
   * random numbers in the interval ]0,1[. To help the phase space
   * generator, the 'dSigHatDR()' should be a smooth function of these
   * numbers, although this is not strictly necessary. Return
   * false if the chosen points failed the kinematical cuts.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(). Uses
   * me().
   */
  virtual CrossSection dSigHatDR() const;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object.
   */
  virtual void setKinematics();
  //@}

  /**
   * Used internally by generateKinematics, after calculating the
   * limits on cos(theta).
   */
  virtual double getCosTheta(double cthmin, double cthmax, const double * r);

  /**
   * Give the option corresponding to the way the scale of an
   * interaction is calculated.
   */
  int scaleChoice() const { return theScaleChoice; }

public:

  /** @name Access cached values in of the last set phase space point. */
  //@{
  /**
   * Return the \f$\hat{t}\f$ of the last set phase space point.
   */
  Energy2 tHat() const { return theLastTHat; }

  /**
   * Return the \f$\hat{u}\f$ of the last set phase space point.
   */
  Energy2 uHat() const { return theLastUHat; }

  /**
   * Return the azimuth angle of the last set phase space point.
   */
  double phi() const { return theLastPhi; }
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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /** @name Set the cached values in of the last set phase space point. */
  //@{
  /**
   * Set the \f$\hat{t}\f$ of the last set phase space point.
   */
  void tHat(Energy2 e2) { theLastTHat = e2; }

  /**
   * Set the \f$\hat{u}\f$ of the last set phase space point.
   */
  void uHat(Energy2 e2) { theLastUHat = e2; }

  /**
   * Set the azimuth angle of the last set phase space point.
   */
  void phi(double phi) { theLastPhi = phi; }
  //@}

protected:

  /**
   * Access to the acutal Switch object used to determine the choice
   * of scale.
   */
  static Switch<ME2to2Base,int> & interfaceScaleChoice();

private:

  /**
   * The option indicating how to calculate the scale of an interaction.
   */
  int theScaleChoice;

  /**
   * The \f$\hat{t}\f$ of the last set phase space point.
   */
  Energy2 theLastTHat;

  /**
   * The \f$\hat{u}\f$ of the last set phase space point.
   */
  Energy2 theLastUHat;

  /**
   * The azimuth angle of the last set phase space point.
   */
  double theLastPhi;


private:

  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<ME2to2Base> initME2to2Base;

  /**
   *  Private and non-existent assignment operator.
   */
  ME2to2Base & operator=(const ME2to2Base &) = delete;

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of ME2to2Base.
 */
template <>
struct BaseClassTrait<ME2to2Base,1>: public ClassTraitsType {
  /** Typedef of the base class of ME2to2Base. */
  typedef MEBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * ME2to2Base class.
 */
template <>
struct ClassTraits<ME2to2Base>: public ClassTraitsBase<ME2to2Base> {
  /** Return the class name. */
  static string className() { return "ThePEG::ME2to2Base"; }
};

/** @endcond */

}

#endif /* ThePEG_ME2to2Base_H */
