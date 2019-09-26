// -*- C++ -*-
//
// MultiCutBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_MultiCutBase_H
#define THEPEG_MultiCutBase_H
//
// This is the declaration of the MultiCutBase class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "MultiCutBase.fh"
#include "Cuts.fh"

namespace ThePEG {

/**
 * This class corresponds to a kinematical cut to be made on a set of
 * outgoing particles from a hard sub-process.
 *
 * There are three virtual functions to be overridden by concrete
 * sub-classes. minS() and maxS() should return the minimum and
 * maximum invariant mass of of a set of particle types. In addition
 * the passCut() function should return true if a set of particle
 * with a given types and given momenta will pass the cuts.
 *
 * @see \ref MultiCutBaseInterfaces "The interfaces"
 * defined for MultiCutBase.
 */
class MultiCutBase: public Interfaced {

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return the minimum allowed value of the squared invariant mass of
   * a set of outgoing partons of the given types. Typically used to
   * cut off the tails of the mass of a resonance for efficiency.
   */
  virtual Energy2 minS(const tcPDVector & pv) const;

  /**
   * Return the maximum allowed value of the squared invariant mass of
   * a set of outgoing partons of the given types. Typically used to
   * cut off the tails of the mass of a resonance for efficiency.
   */
  virtual Energy2 maxS(const tcPDVector & pv) const;

  /**
   * Return true if a set of outgoing particles with typea \a ptype
   * and corresponding momenta \a p passes the cuts.
   */
  virtual bool passCuts(tcCutsPtr parent, const tcPDVector & ptype,
			const vector<LorentzMomentum> & p) const;

  /**
   * Return true if the given vector of particles passes the cuts.
   */
  bool passCuts(tcCutsPtr parent, const tcPVector & p) const;
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
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<MultiCutBase> initMultiCutBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MultiCutBase & operator=(const MultiCutBase &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MultiCutBase. */
template <>
struct BaseClassTrait<MultiCutBase,1> {
  /** Typedef of the first base class of MultiCutBase. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MultiCutBase class and the shared object where it is defined. */
template <>
struct ClassTraits<MultiCutBase>
  : public ClassTraitsBase<MultiCutBase> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::MultiCutBase"; }
};

/** @endcond */

}

#endif /* THEPEG_MultiCutBase_H */
