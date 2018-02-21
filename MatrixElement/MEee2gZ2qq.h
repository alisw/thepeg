// -*- C++ -*-
//
// MEee2gZ2qq.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_MEee2gZ2qq_H
#define ThePEG_MEee2gZ2qq_H
// This is the declaration of the MEee2gZ2qq class.

#include "ThePEG/MatrixElement/ME2to2QCD.h"

namespace ThePEG {

/**
 * The MEee2gZ2qq class implements the
 * \f$e^+e^-\rightarrow\gamma/Z^0\rightarrow q\bar{q}\f$ matrix
 * element. Both the continuum and \f$Z^0\f$ pole term as well as the
 * interference term is included. Although not a strict QCD matrix
 * element the class inherits from ME2to2QCD, mainly to inherit the
 * parameter for the number of active quark flavours.
 *
 * @see \ref MEee2gZ2qqInterfaces "The interfaces"
 * defined for MEee2gZ2qq.
 * @see ME2to2QCD
 */
class MEee2gZ2qq: public ME2to2QCD {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  MEee2gZ2qq();
  //@}

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix element
   * is given. Returns .
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \f$\alpha_{EM}\f$ in which this matrix
   * element is given. Returns 2.
   */
  virtual unsigned int orderInAlphaEW() const;

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

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;
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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

protected:

  /**
   * Constants for the different terms set from the StandardModel in
   * the init() function.
   */
  vector<double> coefs;

  /**
   * The squared mass of the Z0.
   */
  Energy2 mZ2;

  /**
   * The squared width of the Z0.
   */
  Energy2 GZ2;

  /**
   * The last continuum term to be used to select primary diagram.
   */
  mutable double lastCont;

  /**
   * The last Breit-Wigner term to be used to select primary diagram.
   */
  mutable double lastBW;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<MEee2gZ2qq> initMEee2gZ2qq;

  /**
   *  Private and non-existent assignment operator.
   */
  MEee2gZ2qq & operator=(const MEee2gZ2qq &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEee2gZ2qq. */
template <>
struct BaseClassTrait<MEee2gZ2qq,1>: public ClassTraitsType {
  /** Typedef of the first base class of MEee2gZ2qq. */
  typedef ME2to2QCD NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  MEee2gZ2qq class and the shared object where it is defined. */
template <>
struct ClassTraits<MEee2gZ2qq>: public ClassTraitsBase<MEee2gZ2qq> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::MEee2gZ2qq"; }
  /** Return the name of the shared library be loaded to get
   *  access to the MEee2gZ2qq class and every other class it uses
   *  (except the base class). */
  static string library() { return "MEee2gZ2qq.so"; }
};

/** @endcond */

}

#endif /* ThePEG_MEee2gZ2qq_H */
