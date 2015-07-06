// -*- C++ -*-
//
// TwoCutBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_TwoCutBase_H
#define THEPEG_TwoCutBase_H
//
// This is the declaration of the TwoCutBase class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "TwoCutBase.fh"
#include "Cuts.fh"

namespace ThePEG {

/**
 * This class corresponds to a kinematical cut to be made on a pair of
 * particles in a hard sub-process.
 *
 * There are six main virtual functions to be overridden by concrete
 * sub-classes. minsSij(), minTij(), minDeltaR(), minKTClus() and
 * minDurham() returns the minimum allowed values of pre defined
 * kinematical variable. In addition the passCut() function should
 * return true if a pair of particle with a given types and given
 * momenta will pass the cuts.
 *
 * @see \ref TwoCutBaseInterfaces "The interfaces" defined for
 * TwoCutBase.
 */
class TwoCutBase: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  TwoCutBase() {}

  /**
   * The destructor.
   */
  virtual ~TwoCutBase();
  //@}

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return the minimum allowed squared invariant mass of two outgoing
   * partons of type \a pi and \a pj.
   */
  virtual Energy2 minSij(tcPDPtr pi, tcPDPtr pj) const = 0;

  /**
   * Return the minimum allowed value of the negative of the squared
   * invariant mass of an incoming parton of type \a pi and an
   * outgoing parton of type \a po.
   */
  virtual Energy2 minTij(tcPDPtr pi, tcPDPtr po) const = 0;

  /**
   * Return the minimum allowed value of \f$\Delta
   * R_{ij}=\sqrt{\Delta\eta_{ij}^2+\Delta\phi_{ij}^2}\f$ of two
   * outgoing partons of type \a pi and \a pj.
   */
  virtual double minDeltaR(tcPDPtr pi, tcPDPtr pj) const = 0;

  /**
   * Return the minimum allowed value of the longitudinally invariant
   * \f$k_\perp\f$-algorithms distance measure. This is defined as
   * \f$\min(p_{\perp i}, p_{\perp
   * j})\sqrt{\Delta\eta_{ij}^2+\Delta\phi_{ij}^2}\f$ for two outgoing
   * partons, or simply \f$p_{\perp i}\f$ or \f$p_{\perp j}\f$ for a
   * single outgoing parton. Returns 0 if both partons are incoming. A
   * null pointer indicates an incoming parton, hence the type of the
   * incoming parton is irrelevant.
   */
  virtual Energy minKTClus(tcPDPtr pi, tcPDPtr pj) const = 0;

  /**
   * Return the minimum allowed value of the Durham
   * \f$k_\perp\f$-algorithms distance measure. This is defined as
   * \f$2\min(E_j^2, E_j^2)(1-\cos\theta_{ij})/\hat{s}\f$ for two
   * outgoing partons.
   */
  virtual double minDurham(tcPDPtr pi, tcPDPtr pj) const = 0;

  /**
   * Return true if a pair of particles with type \a pitype and \a
   * pjtype and momenta \a pi and \a pj respectively passes the
   * cuts. \a inci and \a inj indicates if the corresponding particles
   * are incoming.
   */
  virtual bool passCuts(tcCutsPtr parent, tcPDPtr pitype, tcPDPtr pjtype,
			LorentzMomentum pi, LorentzMomentum pj,
			bool inci = false, bool incj = false) const;

  /**
   * Return true if the given pair of particles passes the cuts. \a
   * inci and \a inj indicates if the corresponding particles are
   * incoming.
   */
  bool passCuts(tcCutsPtr parent, tcPPtr pi, tcPPtr pj,
		bool inci = false, bool incj = false) const;
  //@}

  /**
   * Describe the currently active cuts in the log file.
   */
  virtual void describe() const;

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractNoPIOClassDescription<TwoCutBase> initTwoCutBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TwoCutBase & operator=(const TwoCutBase &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TwoCutBase. */
template <>
struct BaseClassTrait<TwoCutBase,1> {
  /** Typedef of the first base class of TwoCutBase. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TwoCutBase class and the shared object where it is defined. */
template <>
struct ClassTraits<TwoCutBase>
  : public ClassTraitsBase<TwoCutBase> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::TwoCutBase"; }
};

/** @endcond */

}

#endif /* THEPEG_TwoCutBase_H */
