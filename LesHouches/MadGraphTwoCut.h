// -*- C++ -*-
//
// MadGraphTwoCut.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_MadGraphTwoCut_H
#define THEPEG_MadGraphTwoCut_H
//
// This is the declaration of the MadGraphTwoCut class.
//

#include "ThePEG/Cuts/TwoCutBase.h"

namespace ThePEG {

/**
 * Objects of the MadGraphTwoCut class can be created automatically by
 * the MadGraphReader class when scanning event files for information
 * about cuts. It is also possible to create objects by hand and use
 * it as any other OneCutBase object.
 *
 * @see \ref MadGraphTwoCutInterfaces "The interfaces"
 * defined for MadGraphTwoCut.
 */
class MadGraphTwoCut: public TwoCutBase {

public:

  /**
   * Enumerate the different kinds of cuts made by MadGraph.
   */
  enum class Cut {
    INVMASS, /**< The minimum invariant mass of two particles. */
    DELTAR   /**< The minimum pseudo-rapidity--azimuth-angle distance
                  between two particles. */
  };

  /**
   * Enumerate the types of particles the cut is made on.
   */
  enum class P {
    JET, /**< Coloured particles (jets). */
    LEP, /**< Leptons. */
    PHO, /**< Photons. */
    BOT,  /**< Bottom quarks. */
    NOT  /**< Other types not cut on. */
  };

  /**
   * Enumerate the types of particles pairs the cut is made on.
   */
  enum class PP {
    JETJET, /**< The cut applies only to pairs of coloured particles (jets). */
    LEPLEP, /**< The cut applies only to lepton pairs (in case of INVMASS
                 lepton--anti-lepton pairs of same flavour). */
    PHOPHO, /**< The cut applies only to pairs photons. */
    BOTBOT, /**< The cut applies only to pairs of bottom quarks. */
    BOTJET, /**< The cut applies only to bottom quarks paired with another
                 coloured particle (jet). */
    PHOJET, /**< The cut applies only to a photon paired with a coloured
                 particle (jet). */
    JETLEP, /**< The cut applies only to a coloured particle (jet) paired
                 with a lepton. */
    PHOBOT, /**< The cut applies only to a photon paired with a bottom quark. */
    BOTLEP, /**< The cut applies only to bottom quarks paired with a lepton. */
    PHOLEP  /**< The cut applies only to a photon paired with a lepton. */
  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MadGraphTwoCut()
    : cutType(Cut::DELTAR), pairType(PP::JETJET), theCut(0.0) {}

  /**
   * The constructor used by the MadGraphReader.
   * @param t is the type of the cut.
   * @param p is the type of particles the cut is applied to.
   * @param c is the value of the cut (in units of GeV where applicable).
   */
  MadGraphTwoCut(Cut t, PP p, double c)
    : cutType(t), pairType(p), theCut(c) {}

  //@}

public:

  /** @name Virtual functions mandated by the base class. */
  //@{
  /**
   * Return the minimum allowed squared invariant mass of two outgoing
   * partons of type \a pi and \a pj.
   */
  virtual Energy2 minSij(tcPDPtr pi, tcPDPtr pj) const;

  /**
   * Return the minimum allowed value of the negative of the squared
   * invariant mass of an incoming parton of type \a pi and an
   * outgoing parton of type \a po.
   */
  virtual Energy2 minTij(tcPDPtr pi, tcPDPtr po) const;

  /**
   * Return the minimum allowed value of \f$\Delta
   * R_{ij}=\sqrt{\Delta\eta_{ij}^2+\Delta\phi_{ij}^2}\f$ of two
   * outgoing partons of type \a pi and \a pj.
   */
  virtual double minDeltaR(tcPDPtr pi, tcPDPtr pj) const;

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
  virtual Energy minKTClus(tcPDPtr pi, tcPDPtr pj) const;

  /**
   * Return the minimum allowed value of the Durham
   * \f$k_\perp\f$-algorithms distance measure. This is defined as
   * \f$2\min(E_j^2, E_j^2)(1-\cos\theta_{ij})/\hat{s}\f$ for two
   * outgoing partons.
   */
  virtual double minDurham(tcPDPtr pi, tcPDPtr pj) const;

  /**
   * Return true if a pair of particles with type \a pitype and \a
   * pjtype and momenta \a pi and \a pj respectively passes the
   * cuts. \a inci and \a inj indicates if the corresponding particles
   * are incoming.
   */
  virtual bool passCuts(tcCutsPtr parent, tcPDPtr pitype, tcPDPtr pjtype,
			LorentzMomentum pi, LorentzMomentum pj,
			bool inci = false, bool incj = false) const;
  //@}

protected:

  /**
   * Returns true if cut should be applied to pair of particles of
   * type \a pi and \a pj.
   */
  bool checkType(tcPDPtr pi, tcPDPtr pj) const;

  /**
   * Get the type of particle \a p.
   */
  P getType(tcPDPtr p) const;

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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The type of this cut.
   */
  Cut cutType;

  /**
   * The type of particle pairs this cut applies to.
   */
  PP pairType;

  /**
   * The value of the cut to be applied.
   */
  double theCut;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MadGraphTwoCut> initMadGraphTwoCut;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MadGraphTwoCut & operator=(const MadGraphTwoCut &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MadGraphTwoCut. */
template <>
struct BaseClassTrait<MadGraphTwoCut,1> {
  /** Typedef of the first base class of MadGraphTwoCut. */
  typedef TwoCutBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MadGraphTwoCut class and the shared object where it is defined. */
template <>
struct ClassTraits<MadGraphTwoCut>
  : public ClassTraitsBase<MadGraphTwoCut> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::MadGraphTwoCut"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the MadGraphTwoCut class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "MadGraphReader.so"; }
};

/** @endcond */

}

#endif /* THEPEG_MadGraphTwoCut_H */
