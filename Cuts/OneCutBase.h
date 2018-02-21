// -*- C++ -*-
//
// OneCutBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_OneCutBase_H
#define THEPEG_OneCutBase_H
//
// This is the declaration of the OneCutBase class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "OneCutBase.fh"
#include "Cuts.fh"

namespace ThePEG {

/**
 * This class corresponds to a kinematical cut to be made on a single
 * outgoing parton from a hard sub-process.
 *
 * There are four main virtual functions which must be overridden by
 * concrete sub-classes. minKT() should return the minimum allowed
 * transverse momentum of a given type, while minEta() and maxEta()
 * should return the minimum and maximum allowed pseudo-rapidity for a
 * particle of a given type as measured in the lab-system. Note that
 * when applied in the rest frame of a hard sub-process, the
 * transformation from the lab frame is assumed to be a simple boost
 * along the z-axis. In addition the passCut() function should return
 * true if a particle with a given type and given momentum will pass
 * the cuts.
 *
 * @see \ref OneCutBaseInterfaces "The interfaces"
 * defined for OneCutBase.
 */
class OneCutBase: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  OneCutBase() {}

  /**
   * The destructor.
   */
  virtual ~OneCutBase();
  //@}

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return the minimum allowed value of the transverse momentum of an
   * outgoing parton.
   */
  virtual Energy minKT(tcPDPtr p) const = 0;

  /**
   * Return the minimum allowed pseudo-rapidity of an outgoing parton
   * of the given type. The pseudo-rapidity is measured in the lab
   * system.
   */
  virtual double minEta(tcPDPtr p) const = 0;

  /**
   * Return the maximum allowed pseudo-rapidity of an outgoing parton
   * of the given type. The pseudo-rapidity is measured in the lab
   * system.
   */
  virtual double maxEta(tcPDPtr p) const = 0;

  /**
   * Return the minimum allowed rapidity of an outgoing parton
   * of the given type. The rapidity is measured in the lab
   * system.
   */
  virtual double minRapidityMax(tcPDPtr p) const;

  /**
   * Return the maximum allowed rapidity of an outgoing parton
   * of the given type. The rapidity is measured in the lab
   * system.
   */
  virtual double maxRapidityMin(tcPDPtr p) const;

  /**
   * Return the minimum allowed value of the transverse momentum of
   * the outgoing parton with the lagrest transverse momentum. This
   * version simply returns minKt().
   */
  virtual Energy minMaxKT(tcPDPtr p) const;

  /**
   * Return the minimum allowed pseudo-rapidity of the outgoing parton
   * of the given type with the maximum pseudo-rapidity. The
   * pseudo-rapidity is measured in the lab system. This version
   * simply returns minEta().
   */
  virtual double minMaxEta(tcPDPtr p) const;

  /**
   * Return the maximum allowed pseudo-rapidity of the outgoing parton
   * of the given type with the minimum pseudo-rapidity.. The
   * pseudo-rapidity is measured in the lab system. This version
   * simply returns maxEta().
   */
  virtual double maxMinEta(tcPDPtr p) const;

  /**
   * Return true if a particle with type \a ptype and momentum \a p
   * passes the cuts. The \a parent contains information about the
   * kinematics of the hard sub-process.
   */
  virtual bool passCuts(tcCutsPtr parent,
			tcPDPtr ptype, LorentzMomentum p) const;

  /**
   * Return true if the given particle passes the cuts. The \a parent
   * contains information about the kinematics of the hard
   * sub-process.
   */
  bool passCuts(tcCutsPtr parent, tcPPtr p) const;
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
   * Indicates that this is a concrete class with persistent data.
   */
  static AbstractNoPIOClassDescription<OneCutBase> initOneCutBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  OneCutBase & operator=(const OneCutBase &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of OneCutBase. */
template <>
struct BaseClassTrait<OneCutBase,1> {
  /** Typedef of the first base class of OneCutBase. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the OneCutBase class and the shared object where it is defined. */
template <>
struct ClassTraits<OneCutBase>
  : public ClassTraitsBase<OneCutBase> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::OneCutBase"; }
};

/** @endcond */

}

#endif /* THEPEG_OneCutBase_H */
