// -*- C++ -*-
//
// GRVBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_GRVBase_H
#define ThePEG_GRVBase_H
// This is the declaration of the GRVBase class.

#include "ThePEG/PDF/PDFBase.h"

namespace ThePEG {

/**
 * GRVBase inherits from PDFBase and is used as a base class for all
 * GRV parton densities.
 *
 * @see \ref GRVBaseInterfaces "The interfaces"
 * defined for GRVBase.
 */
class GRVBase: public PDFBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  GRVBase();

  /**
   * Destructor.
   */
  virtual ~GRVBase();
  //@}

public:

  /** @name Virtual functions required by the PDFBase class. */
  //@{
  /**
   * Return true if this PDF can handle the extraction of parton from the
   * given particle, ie. if the particle is a proton or neutron.
   */
  virtual bool canHandleParticle(tcPDPtr particle) const;

  /**
   * Return the parton types which are described by these parton
   * densities.
   */
  virtual cPDVector partons(tcPDPtr p) const;

  /**
   * Return the value of the density of parton at the given a scale
   * and log fractional momentum l (the optional virtuality of the
   * incoming particle is not used).
   */
  virtual double xfl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale) const;

  /**
   * Return the valaens partof the density of parton at the given a
   * scale and log fractional momentum l (the optional virtuality of
   * the incoming particle is not used).
   */
  virtual double xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale) const;
  //@}

public:

  /** @name Access derived kinematical quantities. */
  //@{
  /**
   * Return last selected
   * \f$S\f$. \f$S=\log(\log(Q^2/\mu^2)/\log(Q^2/\Lambda_{QCD}^2))\f$
   */
  double S() const { return theS; }

  /**
   * Return last selected
   * \f$S^2\f$. \f$S=\log(\log(Q^2/\mu^2)/\log(Q^2/\Lambda_{QCD}^2))\f$
   */
  double S2() const { return theS2; }

  /**
   * Return last selected
   * \f$S^3\f$. \f$S=\log(\log(Q^2/\mu^2)/\log(Q^2/\Lambda_{QCD}^2))\f$
   */
  double S3() const { return theS3; }

  /**
   * Return last selected
   * \f$\sqrt{S}\f$. \f$S=\log(\log(Q^2/\mu^2)/\log(Q^2/\Lambda_{QCD}^2))\f$
   */
  double rootS() const { return theRootS; }

  /**
   * Return last selected momentum fraction, \f$x\f$.
   */
  double x() const { return thex; }

  /**
   * Return last selected logarithmic momentum fraction
   * \f$l=\log(1/x)\f$.
   */
  double lx() const { return theLx; }

  /**
   * Return one minus the last selected momentum fraction, eps\f$=1-x\f$.
   */
  double eps() const { return theEps; }

  /**
   * Return the square root of the last selected momentum fraction,
   * \f$x\f$.
   */
  double rootx() const { return theRootx; }

  //@}

protected:

  /**
   * Setup the \a l\f$=\log{1/x}\f$ and \a scale \f$Q^2\f$ to be used
   * in the following call to uv(), dv)=, etc.
   */
  virtual void setup(double l, Energy2 scale) const = 0;

  /**
   * Setup the \a l\f$=\log{1/x}\f$ and \a scale \f$Q^2\f$ to be used
   * in the following call to uv(), dv)=, etc.
   */
  void setup(double l, Energy2 scale, Energy2 mu2, Energy2 lam2) const;

  /**
   * The form of the valens density functions.
   */
  double valens(double N, double ak, double bk,
		double a, double b, double c, double d) const;

  /**
   * The form of the light sea and gluon density
   * functions.
   */
  double lightsea(double al, double be, double ak, double bk, double a,
		  double b, double c, double d, double e, double es) const;

  /**
   * The form of the heavy sea density functions.
   */
  double heavysea(double sth, double al, double be, double ak, double ag,
		  double b, double d, double e, double es) const;

  /**
   * Return the value of the u valens density for the values previously given
   * by setup().
   */
  virtual double uv() const = 0;

  /**
   * Return the value of the d valens density for the values previously given
   * by setup().
   */
  virtual double dv() const = 0;

  /**
   * Return the value of the difference between the u and d sea
   * densities for the values previously given by setup().
   */
  virtual double del() const = 0;

  /**
   * Return the value of the average u and d sea densities for the
   * values previously given by setup().
   */
  virtual double udb() const = 0;

  /**
   * Return the value of the s density for the values previously given by
   * setup().
   */
  virtual double sb() const = 0;

  /**
   * Return the value of the c density for the values previously given by
   * setup().
   */
  virtual double cb() const = 0;

  /**
   * Return the value of the b density for the values previously given by
   * setup().
   */
  virtual double bb() const = 0;

  /**
   * Return the value of the gluon densities for the values previously
   * given by setup().
   */
  virtual double gl() const = 0;

  /**
   * fuv() returns the saved values from the quv() functions if
   * present. Otherwise uv() is called, saved and returned.
   */
  double fuv() const { return uvSave >= 0.0? uvSave: ( uvSave = uv() ); }

  /**
   * fdv() returns the saved values from the dv() functions if
   * present. Otherwise dv() is called, saved and returned.
   */
  double fdv() const { return dvSave >= 0.0? dvSave: ( dvSave = dv() ); }

  /**
   * fdel() returns the saved values from the del() functions if
   * present. Otherwise del() is called, saved and returned.
   */
  double fdel() const { return delSave >= 0.0? delSave: ( delSave = del() ); }

  /**
   * fudb() returns the saved values from the udb() functions if
   * present. Otherwise udb() is called, saved and returned.
   */
  double fudb() const { return udbSave >= 0.0? udbSave: ( udbSave = udb() ); }

  /**
   * fsb() returns the saved values from the sb() functions if
   * present. Otherwise sb() is called, saved and returned.
   */
  double fsb() const { return sbSave >= 0.0? sbSave: ( sbSave = sb() ); }

  /**
   * fcb() returns the saved values from the cb() functions if
   * present. Otherwise cb() is called, saved and returned.
   */
  double fcb() const { return cbSave >= 0.0? cbSave: ( cbSave = cb() ); }

  /**
   * fbb() returns the saved values from the bb() functions if
   * present. Otherwise bb() is called, saved and returned.
   */
  double fbb() const { return bbSave >= 0.0? bbSave: ( bbSave = bb() ); }

  /**
   * fgl() returns the saved values from the gl() functions if
   * present. Otherwise gl() is called, saved and returned.
   */
  double fgl() const { return glSave >= 0.0? glSave: ( glSave = gl() ); }

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

private:

  /**
   * The last selected logarithmic momentum fraction
   * \f$l=\log(1/x)\f$.
   */
  mutable double theLx;

  /**
   * THe last selected momentum fraction, \f$x\f$.
   */
  mutable double thex;

  /**
   * One minus the last selected momentum fraction, eps\f$=1-x\f$.
   */
  mutable double theEps;

  /**
   * The square root of the last selected momentum fraction, \f$x\f$.
   */
  mutable double theRootx;

  /**
   * The last selected scale.
   */
  mutable Energy2 Q2;

  /**
   * The last used \f$\Lambda_{QCD}^2\f$.
   */
  mutable Energy2 theLam2;

  /**
   * The last used \f$\mu^2\f$.
   */
  mutable Energy2 theMu2;

  /**
   * The last selected
   * \f$S\f$. \f$S=\log(\log(Q^2/\mu^2)/\log(Q^2/\Lambda_{QCD}^2))\f$
   */
  mutable double theS;

  /**
   * Return last selected \f$S^2\f$.
   */
  mutable double theS2;

  /**
   * Return last selected \f$S^3\f$.
   */
  mutable double theS3;

  /**
   * Return last selected \f$\sqrt{S}\f$.
   */
  mutable double theRootS;

  /**
   * Saved values from the different functions.
   */
  mutable double uvSave;

  /**
   * Saved values from the different functions.
   */
  mutable double dvSave;

  /**
   * Saved values from the different functions.
   */
  mutable double delSave;

  /**
   * Saved values from the different functions.
   */
  mutable double udbSave;

  /**
   * Saved values from the different functions.
   */
  mutable double sbSave;

  /**
   * Saved values from the different functions.
   */
  mutable double cbSave;

  /**
   * Saved values from the different functions.
   */
  mutable double bbSave;

  /**
   * Saved values from the different functions.
   */
  mutable double glSave;

private:

  /**
   *  Private and non-existent assignment operator.
   */
  GRVBase & operator=(const GRVBase &) = delete;

};

}

#endif /* ThePEG_GRVBase_H */
