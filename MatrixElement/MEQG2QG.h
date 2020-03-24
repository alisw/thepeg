// -*- C++ -*-
//
// MEQG2QG.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_MEQG2QG_H
#define ThePEG_MEQG2QG_H
// This is the declaration of the MEQG2QG class.

#include "ThePEG/MatrixElement/ME2to2QCD.h"

namespace ThePEG {

/**
 * MEQG2QG inherits from ME2to2QCD and implements the
 * standard \f$qg\rightarrow qg\f$ matrix element.
 *
 * @see \ref MEQG2QGInterfaces "The interfaces"
 * defined for MEQG2QG.
 * @see ME2to2QCD
 */
class MEQG2QG: public ME2to2QCD {

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;
  //@}

protected:

  /** @name Internal functions returning the matrix element squared
   *  for different colour configurations. */
  //@{
  /**
   * Return the matrix element squared (without common pre-factors)
   * for the specific colour configuration.
   */
  double colA1() const 
  { 
    return ( interference()? 2.25: 2.0 )*sqr(uHat()/tHat()); 
  }

  /**
   * Return the matrix element squared (without common pre-factors)
   * for the specific colour configuration.
   */
  double colB1() const 
  { 
    return ( interference()? 2.25: 2.0 )*sqr(sHat()/tHat());
  }

  /**
   * Return the matrix element squared (without common pre-factors)
   * for the specific colour configuration.
   */
  double colA2() const { return  -uHat()/sHat(); }

  /**
   * Return the matrix element squared (without common pre-factors)
   * for the specific colour configuration.
   */
  double colB2() const { return -sHat()/uHat(); }
  //@}

public:

  /**
   * Standard Init function used to initialize the interfaces.
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
  static NoPIOClassDescription<MEQG2QG> initMEQG2QG;

  /**
   *  Private and non-existent assignment operator.
   */
  MEQG2QG & operator=(const MEQG2QG &) = delete;

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEQG2QG. */
template <>
struct BaseClassTrait<MEQG2QG,1>: public ClassTraitsType {
  /** Typedef of the first base class of MEQG2QG. */
  typedef ME2to2QCD NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEQG2QG class and the shared object where it is defined. */
template <>
struct ClassTraits<MEQG2QG>: public ClassTraitsBase<MEQG2QG> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::MEQG2QG"; }
  /** Return the name of the shared library be loaded to get
   *  access to the MEQG2QG class and every other class it uses
   *  (except the base class). */
  static string library() { return "MEQCD.so"; }
};

/** @endcond */

}

#endif /* ThePEG_MEQG2QG_H */
