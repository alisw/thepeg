// -*- C++ -*-
//
// V2LeptonsCut.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_V2LeptonsCut_H
#define THEPEG_V2LeptonsCut_H
//
// This is the declaration of the V2LeptonsCut class.
//

#include "ThePEG/Cuts/MultiCutBase.h"

namespace ThePEG {

/**
 * This class inherits from MultiCutBase and describes cuts on the
 * invariant mass of two final state leptons corresponding to the
 * decay of a vector boson.  It can be used when generating matrix
 * elements to avoid the long tails of the resonance.
 *
 * @see \ref V2LeptonsCutInterfaces "The interfaces"
 * defined for V2LeptonsCut.
 */
class V2LeptonsCut: public MultiCutBase {

  /**
   * Enumeration of the different families.
   */
  enum Family {
    electron = 1, /**< Lepton Family. */
    muon = 2,     /**< Muon Family. */
    tau = 4       /**< Tau Family. */
  };

  /**
   * Enumeration of charge combinations.
   */
  enum CComb {
    posneg = 1, /**< charged lepton anti-lepton pair. */
    negneu = 2, /**< negative lepton anti-neutrino pair. */
    posneu = 4, /**< positive lepton anti-neutrino pair. */
    neuneu = 8  /**< neutrino anti-neutrino pair. */
  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  V2LeptonsCut() : theMinM(70.0*GeV), theMaxM(90.0*GeV), theFamilies(electron|muon),
		   theCComb(negneu|posneu) {}

  /**
   * The destructor.
   */
  virtual ~V2LeptonsCut();
  //@}

public:

  /** @name Overridden virtual functions defined in the base class. */
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
  //@}

  /**
   * Describe the currently active cuts in the log file.
   */
  virtual void describe() const;

protected:

  /**
   * Check if the PDG id numbers matches this cut.
   */
  bool checkTypes(long id1, long id2) const;

  /**
   * Check the family of the given PDG id number.
   */
  int family(long id) const;

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
  Energy maxMinM() const;

  /**
   * Helper function used by the interface.
   */
  Energy minMaxM() const;

private:

  /**
   * The minimum invariant mass.
   */
  Energy theMinM;

  /**
   * The maximum invariant mass.
   */
  Energy theMaxM;

  /**
   * Integer corresponding to the lepton families to match.
   */
  int theFamilies;

  /**
   * Integer corresponding to the charge combination to match.
   */
  int theCComb;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<V2LeptonsCut> initV2LeptonsCut;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  V2LeptonsCut & operator=(const V2LeptonsCut &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of V2LeptonsCut. */
template <>
struct BaseClassTrait<V2LeptonsCut,1> {
  /** Typedef of the first base class of V2LeptonsCut. */
  typedef MultiCutBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the V2LeptonsCut class and the shared object where it is defined. */
template <>
struct ClassTraits<V2LeptonsCut>
  : public ClassTraitsBase<V2LeptonsCut> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::V2LeptonsCut"; }
  /** Return the name of the shared library be loaded to get
   *  access to the V2LeptonsCut class and every other class it uses
   *  (except the base class). */
  static string library() { return "V2LeptonsCut.so"; }
};

/** @endcond */

}

#endif /* THEPEG_V2LeptonsCut_H */
