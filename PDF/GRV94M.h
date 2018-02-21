// -*- C++ -*-
//
// GRV94M.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_GRV94M_H
#define ThePEG_GRV94M_H
// This is the declaration of the GRV94M class.

#include "ThePEG/PDF/GRVBase.h"

namespace ThePEG {

/**
 * GRV94M iherits from PDFBase via the GRVBase class and implements
 * the GRV94M parton densities for (anti) protons ad neutrons.
 *
 * @see \ref GRV94MInterfaces "The interfaces"
 * defined for GRV94M.
 */
class GRV94M: public GRVBase {

  /**
   * Return the cutoff scale.
   */
  Energy2 mu2() const { return 0.34*GeV2; }

  /**
   * Return the square of \f$\Lambda_{QCD}\f$ used.
   */
  Energy2 lam2() const { return sqr(0.248*GeV); }

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
  static NoPIOClassDescription<GRV94M> initGRV94M;

  /**
   *  Private and non-existent assignment operator.
   */
  GRV94M & operator=(const GRV94M &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GRV94M. */
template <>
struct BaseClassTrait<GRV94M,1>: public ClassTraitsType {
  /** Typedef of the first base class of GRV94M. */
  typedef GRVBase NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  GRV94M class and the shared object where it is defined. */
template <>
struct ClassTraits<GRV94M>: public ClassTraitsBase<GRV94M> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::GRV94M"; }
  /** Return the name of the shared library be loaded to get access to
   *  the GRV94M class and every other class it uses (except
   *  the base class). */
  static string library() { return "GRV94M.so"; }
};

/** @endcond */

}

#endif /* ThePEG_GRV94M_H */
