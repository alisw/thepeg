// -*- C++ -*-
//
// LWHFactory.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_LWHFactory_H
#define THEPEG_LWHFactory_H
//
// This is the declaration of the LWHFactory class.
//

#include "ThePEG/Analysis/FactoryBase.h"

namespace ThePEG {

/**
 * Here is the documentation of the LWHFactory class. It inherits from
 * the abstract FactoryBase class and implements the Light-Weight
 * Histogram package, LWH. This implements the most rudimentary
 * histogramming facilities according to the <a
 * href="http://aida.freehep.org">AIDA</a> interface
 * specifications. Currently the only thing that is supported is
 * simple, equally binned, one dimensional histograms. It is mainly
 * intended to be used in applications where one needs to fill simple
 * histograms and output them. With LWH it is then possible to do this
 * without the overhead of a full AIDA implementation, but still
 * having the option to use a full implementation later on with
 * minimal changes.
 *
 * @see \ref LWHFactoryInterfaces "The interfaces"
 * defined for LWHFactory.
 */
class LWHFactory: public FactoryBase {

public:

  /** @name Manipulate histograms */
  //@{
  /**
   * Rescale the given \a histogram so that the integral over the bins
   * will give the correct integrated cross section for the observable
   * in the given \a unit.
   */
  virtual void
  normalizeToXSec(tH1DPtr histogram, CrossSection unit = picobarn) const;

  /**
   * Rescale the given \a histogram so that the integral over the bins
   * gives the fraction of the total cross section generated which is
   * contained in the bins.
   */
  virtual void normalizeToXSecFraction(tH1DPtr histogram) const;

  /**
   * Rescale the given \a histogram so that the integral over the bins
   * gives one.
   */
  virtual void normalizeToUnity(tH1DPtr histogram) const;

  /**
   * Rescale the given \a histogram so that the integral over the bins
   * will give the correct integrated cross section for the observable
   * in the given \a unit.
   */
  virtual void
  normalizeToXSec(tH2DPtr histogram, CrossSection unit = picobarn) const;

  /**
   * Rescale the given \a histogram so that the integral over the bins
   * gives the fraction of the total cross section generated which is
   * contained in the bins.
   */
  virtual void normalizeToXSecFraction(tH2DPtr histogram) const;

  /**
   * Rescale the given \a histogram so that the integral over the bins
   * gives one.
   */
  virtual void normalizeToUnity(tH2DPtr histogram) const;
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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const { return new_ptr(*this); }
  //@}



protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<LWHFactory> initLWHFactory;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LWHFactory & operator=(const LWHFactory &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LWHFactory. */
template <>
struct BaseClassTrait<LWHFactory,1> {
  /** Typedef of the first base class of LWHFactory. */
  typedef FactoryBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LWHFactory class and the shared object where it is defined. */
template <>
struct ClassTraits<LWHFactory>
  : public ClassTraitsBase<LWHFactory> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::LWHFactory"; }

  /** Return the name of the shared library be loaded to get access to
   *  the LWHFactory class and every other class it uses
   *  (except the base class). */
  static string library() { return "LWHFactory.so"; }
};

/** @endcond */

}

#endif /* THEPEG_LWHFactory_H */
