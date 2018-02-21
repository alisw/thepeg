// -*- C++ -*-
//
// SimpleKTCut.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_SimpleKTCut_H
#define THEPEG_SimpleKTCut_H
//
// This is the declaration of the SimpleKTCut class.
//

#include "ThePEG/Cuts/OneCutBase.h"

namespace ThePEG {

/**
 * This is a very simple concrete sub-class of OneCutbase simply
 * requiring a minimum transverse momentum of any outgoing
 * particle. It is also possible to require a minimum and maximum
 * pseudorapidity. Optionally the restrictions only apply to particles
 * matching a specific matcher object.
 *
 * @see \ref SimpleKTCutInterfaces "The interfaces"
 * defined for SimpleKTCut.
 */
class SimpleKTCut: public OneCutBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SimpleKTCut(Energy minKT=10*GeV) 
    : theMinKT(minKT), theMaxKT(Constants::MaxEnergy),
      theMinEta(-Constants::MaxRapidity),
      theMaxEta(Constants::MaxRapidity) {}
  
  /**
   * The destructor.
   */
  virtual ~SimpleKTCut();
  //@}

public:

  /** @name Overwritten virtual functions defined in the base class. */
  //@{
  /**
   * Return the minimum allowed value of the transverse momentum of an
   * outgoing parton.
   */
  virtual Energy minKT(tcPDPtr p) const;

  /**
   * Return the minimum allowed pseudo-rapidity of an outgoing parton
   * of the given type. The pseudo-rapidity is measured in the lab
   * system.
   */
  virtual double minEta(tcPDPtr p) const;

  /**
   * Return the maximum allowed pseudo-rapidity of an outgoing parton
   * of the given type. The pseudo-rapidity is measured in the lab
   * system.
   */
  virtual double maxEta(tcPDPtr p) const;

  /**
   * Return true if a particle with type \a ptype and momentum \a p
   * passes the cuts. The \a parent contains information about the
   * kinematics of the hard sub-process.
   */
  virtual bool passCuts(tcCutsPtr parent,
			tcPDPtr ptype, LorentzMomentum p) const;
  //@}

  /**
   * Describe the currently active cuts in the log file.
   */
  virtual void describe() const;

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
   * Helper function used by the interface.
   */
  Energy maxKTMin() const;

  /**
   * Helper function used by the interface.
   */
  Energy minKTMax() const;

  /**
   * Helper function used by the interface.
   */
  double maxEtaMin() const;

  /**
   * Helper function used by the interface.
   */
  double minEtaMax() const;

private:

  /**
   * The minimum allowed value of the transverse momentum of an
   * outgoing parton.
   */
  Energy theMinKT;

  /**
   * The maximum allowed value of the transverse momentum of an
   * outgoing parton.
   */
  Energy theMaxKT;

  /**
   * The minimum allowed pseudo-rapidity of an outgoing parton. The
   * pseudo-rapidity is measured in the lab system.
   */
  double theMinEta;

  /**
   * The maximum allowed pseudo-rapidity of an outgoing parton. The
   * pseudo-rapidity is measured in the lab system.
   */
  double theMaxEta;

  /**
   * If non-null only particles matching this object will be affected
   * by this cut.
   */
  PMPtr theMatcher;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SimpleKTCut> initSimpleKTCut;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SimpleKTCut & operator=(const SimpleKTCut &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SimpleKTCut. */
template <>
struct BaseClassTrait<SimpleKTCut,1> {
  /** Typedef of the first base class of SimpleKTCut. */
  typedef OneCutBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SimpleKTCut class and the shared object where it is defined. */
template <>
struct ClassTraits<SimpleKTCut>
  : public ClassTraitsBase<SimpleKTCut> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::SimpleKTCut"; }
  /** Return the name of the shared library be loaded to get
   *  access to the SimpleKTCut class and every other class it uses
   *  (except the base class). */
  static string library() { return "SimpleKTCut.so"; }
};

/** @endcond */

}

#endif /* THEPEG_SimpleKTCut_H */
