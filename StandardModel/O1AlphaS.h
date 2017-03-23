// -*- C++ -*-
//
// O1AlphaS.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_O1AlphaS_H
#define ThePEG_O1AlphaS_H
// This is the declaration of the O1AlphaS class.

#include "AlphaSBase.h"

namespace ThePEG {

/**
 * O1AlphaS inherits from AlphaSBase and implements the leading order
 * running QCD coupling. The value is determined by a
 * \f$\Lambda_{QCD}\f$ parameter at a given number of
 * flavours. Optionally the coupling can be frozen under some minimum
 * scale to avoid divergencies or negative couplings.
 *
 * @see \ref O1AlphaSInterfaces "The interfaces"
 * defined for O1AlphaS.
 */
class O1AlphaS: public AlphaSBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  O1AlphaS()
    : theLambdaQCD(0.25*GeV), theLambdaFlavour(4), 
      theMaxFlav(6), Q0(ZERO) {}
  //@}

public:

  /** @name Virtual functions mandated by the sub-class. */
  //@{
  /**
   * The \f$\alpha_S\f$. Return the QCD coupling for a given \a scale
   * using the given standard model object \a sm.
   */
  virtual double value(Energy2 scale, const StandardModelBase &) const;

  /**
   * Return the number of loops contributing to
   * the running this coupling.
   */
  virtual unsigned int nloops () const { return 1; }

  /**
   * Return the flavour thresholds used. The returned vector contains
   * (in position <code>i</code>) the scales when the active number of
   * flavours changes from <code>i</code> to <code>i+1</code>.
   */
  virtual vector<Energy2> flavourThresholds() const;

  /**
   * Return the \f$\Lambda_{QCD}\f$ used for different numbers of
   * active flavours.
   */
  virtual vector<Energy> LambdaQCDs() const;
  //@}

  /**
   * Return the maximum number of active flavours.
   */
  int getMaxFlav() const { return theMaxFlav; }

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


private:

  /**
   * The \f$\Lambda_{QCD}\f$ for the number of flavours specified by
   * theLambdaFlavour. Other \f$\Lambda_{QCD}\f$ values for other
   * numbers of active flavours are calculated from
   * flavourThresholds() using a continuity requirement.
   */
  Energy theLambdaQCD;

  /**
   * The number of flavours for which theLambdaQCD is given.
   */
  int theLambdaFlavour;

  /**
   * The maximum number of active flavours.
   */
  int theMaxFlav;

  /**
   * The scale below which \f$\alpha_S\f$ is frozen.
   */
  Energy Q0;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<O1AlphaS> initO1AlphaS;

  /**
   *  Private and non-existent assignment operator.
   */
  O1AlphaS & operator=(const O1AlphaS &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of O1AlphaS. */
template <>
struct BaseClassTrait<O1AlphaS,1>: public ClassTraitsType {
  /** Typedef of the first base class of O1AlphaS. */
  typedef AlphaSBase NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  O1AlphaS class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<O1AlphaS>: public ClassTraitsBase<O1AlphaS> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::O1AlphaS"; }
  /** Return the name of the shared library be loaded to get access to
   *  the O1AlphaS class and every other class it uses
   *  (except the base class). */
  static string library() { return "O1AlphaS.so"; }
};

/** @endcond */

}

#endif /* ThePEG_O1AlphaS_H */
