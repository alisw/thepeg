// -*- C++ -*-
//
// GRV94L.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_GRV94L_H
#define ThePEG_GRV94L_H
// This is the declaration of the GRV94L class.

#include "ThePEG/PDF/GRVBase.h"

namespace ThePEG {

/**
 * GRV94L inherits from PDFBase via the GRVBase class and implements
 * the GRV94L parton densities for (anti) protons and neutrons.
 *
 * @see \ref GRV94LInterfaces "The interfaces"
 * defined for GRV94L.
 */
class GRV94L: public GRVBase {

  /**
   * Return the cutoff scale.
   */
  Energy2 mu2() const { return 0.23*GeV2; }

  /**
   * Return the square of \f$\Lambda_{QCD}\f$ used.
   */
  Energy2 lam2() const { return sqr(0.2322*GeV); }

protected:

  /**
   * Setup the \a l\f$=\log{1/x}\f$ and \a scale \f$Q^2\f$ to be used
   * in the following call to uv(), dv)=, etc.
   */
  virtual void setup(double l, Energy2 scale) const;

  /**
   * Return the value of the u valens density for the values previously given
   * by setup().
   */
  virtual double uv() const;

  /**
   * Return the value of the d valens density for the values previously given
   * by setup().
   */
  virtual double dv() const;

  /**
   * Return the value of the difference between the u and d sea
   * densities for the values previously given by setup().
   */
  virtual double del() const;

  /**
   * Return the value of the average u and d sea densities for the
   * values previously given by setup().
   */
  virtual double udb() const;

  /**
   * Return the value of the s density for the values previously given by
   * setup().
   */
  virtual double sb() const;

  /**
   * Return the value of the c density for the values previously given by
   * setup().
   */
  virtual double cb() const;

  /**
   * Return the value of the b density for the values previously given by
   * setup().
   */
  virtual double bb() const;

  /**
   * Return the value of the gluon densities for the values previously
   * given by setup().
   */
  virtual double gl() const;

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
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<GRV94L> initGRV94L;

  /**
   *  Private and non-existent assignment operator.
   */
  GRV94L & operator=(const GRV94L &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GRV94L. */
template <>
struct BaseClassTrait<GRV94L,1>: public ClassTraitsType {
  /** Typedef of the first base class of GRV94L. */
  typedef GRVBase NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  GRV94L class and the shared object where it is defined. */
template <>
struct ClassTraits<GRV94L>: public ClassTraitsBase<GRV94L> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::GRV94L"; }
  /** Return the name of the shared library be loaded to get access to
   *  the GRV94L class and every other class it uses (except
   *  the base class). */
  static string library() { return "GRV94L.so"; }
};

/** @endcond */

}

#endif /* ThePEG_GRV94L_H */
