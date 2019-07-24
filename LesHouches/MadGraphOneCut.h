// -*- C++ -*-
//
// MadGraphOneCut.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_MadGraphOneCut_H
#define THEPEG_MadGraphOneCut_H
//
// This is the declaration of the MadGraphOneCut class.
//

#include "ThePEG/Cuts/OneCutBase.h"

namespace ThePEG {

/**
 * Objects of the MadGraphOneCut class can be created automatically by
 * the MadGraphReader class when scanning event files for information
 * about cuts. It is also possible to create objects by hand and use
 * it as any other OneCutBase object.
 *
 * @see \ref MadGraphOneCutInterfaces "The interfaces"
 * defined for MadGraphOneCut.
 */
class MadGraphOneCut: public OneCutBase {

public:

  /**
   * Enumerate the different kinds of cuts made by MadGraph.
   */
  enum class Cut {
    PT,  /**< The minimum transverse momentum of a particle. */
    ETA, /**< The maximum (absolute value of) pseudo-rapidity of a particle. */
    XPT  /**< The minimum transverse momentum of the particle with
	      largest transverse momentum. */
  };

  /**
   * Enumerate the types of particles the cut is made on.
   */
  enum class P {
    JET, /**< The cut applies only to coloured particles. */
    LEP, /**< The cut applies only to leptons. */
    PHO, /**< The cut applies only to photons. */
    BOT  /**< The cut applies only to bottom quarks. */
  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MadGraphOneCut() : cutType(Cut::PT), particleType(P::JET), theCut(0.0) {}

  /**
   * The constructor used by the MadGraphReader.
   * @param t is the type of the cut.
   * @param p is the type of particles the cut is applied to.
   * @param c is the value of the cut (in units of GeV where applicable).
   */
  MadGraphOneCut(Cut t, P p, double c)
    : cutType(t), particleType(p), theCut(c) {}
  //@}

public:

  /** @name Virtual functions mandated by the base class. */
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
   * Return the minimum allowed value of the transverse momentum of
   * the outgoing parton with the lagrest transverse momentum. This
   * version simply returns minKt().
   */
  virtual Energy minMaxKT(tcPDPtr p) const;

  /**
   * Return true if a particle with type \a ptype and momentum \a p
   * passes the cuts. The \a parent contains information about the
   * kinematics of the hard sub-process.
   */
  virtual bool passCuts(tcCutsPtr parent,
			tcPDPtr ptype, LorentzMomentum p) const;
  //@}

protected:

  /**
   * Returns true if cut should be applied to a particle of type \a p.
   */
  bool checkType(tcPDPtr p) const;

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
   * The type of particles this cut applies to.
   */
  P particleType;

  /**
   * The value of the cut to be applied.
   */
  double theCut;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MadGraphOneCut> initMadGraphOneCut;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MadGraphOneCut & operator=(const MadGraphOneCut &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MadGraphOneCut. */
template <>
struct BaseClassTrait<MadGraphOneCut,1> {
  /** Typedef of the first base class of MadGraphOneCut. */
  typedef OneCutBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MadGraphOneCut class and the shared object where it is defined. */
template <>
struct ClassTraits<MadGraphOneCut>
  : public ClassTraitsBase<MadGraphOneCut> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::MadGraphOneCut"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the MadGraphOneCut class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "MadGraphReader.so"; }
};

/** @endcond */

}

#endif /* THEPEG_MadGraphOneCut_H */
