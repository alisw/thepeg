// -*- C++ -*-
//
// NJetsCut.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
// Copyright (C) 2009-2017 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_NJetsCut_H
#define THEPEG_NJetsCut_H
//
// This is the declaration of the NJetsCut class.
//

#include "ThePEG/Cuts/MultiCutBase.h"
#include "ThePEG/PDT/MatcherBase.h"

namespace ThePEG {

/**
 * NJetsCut is a simple cut on jet multiplicity.
 *
 * @see \ref NJetsCutInterfaces "The interfaces"
 * defined for NJetsCut.
 */
class NJetsCut: public MultiCutBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  NJetsCut();

  /**
   * The destructor.
   */
  virtual ~NJetsCut();
  //@}

public:

  /** @name Overridden virtual functions defined in the base class. */
  //@{
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
   * A matcher for unresolved partons.
   */
  Ptr<MatcherBase>::ptr unresolvedMatcher;

  /**
   * The minimum number of jets.
   */
  int nJetsMin;

  /**
   * The maximum number of jets.
   */
  int nJetsMax;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<NJetsCut> initNJetsCut;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NJetsCut & operator=(const NJetsCut &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NJetsCut. */
template <>
struct BaseClassTrait<NJetsCut,1> {
  /** Typedef of the first base class of NJetsCut. */
  typedef MultiCutBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NJetsCut class and the shared object where it is defined. */
template <>
struct ClassTraits<NJetsCut>
  : public ClassTraitsBase<NJetsCut> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::NJetsCut"; }
  /** Return the name of the shared library be loaded to get
   *  access to the NJetsCut class and every other class it uses
   *  (except the base class). */
  static string library() { return "JetCuts.so"; }
};

/** @endcond */

}

#endif /* THEPEG_NJetsCut_H */
