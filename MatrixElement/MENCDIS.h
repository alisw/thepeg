// -*- C++ -*-
//
// MENCDIS.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_MENCDIS_H
#define ThePEG_MENCDIS_H
// This is the declaration of the MENCDIS class.

#include "ThePEG/MatrixElement/ME2to2QCD.h"
// #include "MENCDIS.fh"
// #include "MENCDIS.xh"

namespace ThePEG {

/**
 * The MENCDIS class implements the \f$e^\pm q\rightarrow e^\pm q\f$
 * matrix element. Both the gamma and \f$Z^0\f$ terms as well
 * as the interference term is included. Although not a strict QCD
 * matrix element the class inherits from ME2to2QCD, mainly to inherit
 * the parameter for the number of active quark flavours.
 *
 * @see \ref MENCDISInterfaces "The interfaces"
 * defined for MENCDIS.
 * @see ME2to2QCD
 */
class MENCDIS: public ME2to2QCD {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  MENCDIS();

  /**
   * Copy-constructor.
   */
  MENCDIS(const MENCDIS &);

  /**
   * Destructor.
   */
  virtual ~MENCDIS();
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
   * The squared mass of the Z0.
   */
  Energy2 mZ2;

  /**
   * The last gamma term to be used to select primary diagram.
   */
  mutable double lastG;

  /**
   * The last Z0 term to be used to select primary diagram.
   */
  mutable double lastZ;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<MENCDIS> initMENCDIS;

  /**
   *  Private and non-existent assignment operator.
   */
  MENCDIS & operator=(const MENCDIS &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MENCDIS. */
template <>
struct BaseClassTrait<MENCDIS,1>: public ClassTraitsType {
  /** Typedef of the first base class of MENCDIS. */
  typedef ME2to2QCD NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  MENCDIS class and the shared object where it is defined. */
template <>
struct ClassTraits<MENCDIS>: public ClassTraitsBase<MENCDIS> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::MENCDIS"; }
  /** Return the name of the shared library be loaded to get
   *  access to the MENCDIS class and every other class it uses
   *  (except the base class). */
  static string library() { return "MENCDIS.so"; }
};

/** @endcond */

}

#endif /* ThePEG_MENCDIS_H */
