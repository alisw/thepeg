// -*- C++ -*-
//
// JetCuts.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
// Copyright (C) 2009-2017 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_JetCuts_H
#define THEPEG_JetCuts_H
//
// This is the declaration of the JetCuts class.
//

#include "ThePEG/Cuts/MultiCutBase.h"
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/Cuts/JetRegion.h"
#include "ThePEG/Cuts/JetPairRegion.h"
#include "ThePEG/Cuts/MultiJetRegion.h"

namespace ThePEG {

/**
 * JetCuts combines various JetRegion and JetPairRegion objects into a
 * cut object.
 *
 * @see JetRegion
 * @see JetPairRegion
 *
 * @see \ref JetCutsInterfaces "The interfaces"
 * defined for JetCuts.
 */
class JetCuts: public MultiCutBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  JetCuts();

  /**
   * The destructor.
   */
  virtual ~JetCuts();
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

  /**
   * Return the matcher for unresolved partons.
   */
  Ptr<MatcherBase>::tptr unresolvedMatcher() const { return theUnresolvedMatcher; }

  /**
   * Return the jet regions to match.
   */
  const vector<Ptr<JetRegion>::ptr>& jetRegions() const { return theJetRegions; }

  /**
   * Return the jet veto regions to check.
   */
  const vector<Ptr<JetRegion>::ptr>& jetVetoRegions() const { return theJetVetoRegions; }

  /**
   * Return the jet pair regions to match.
   */
  const vector<Ptr<JetPairRegion>::ptr>& jetPairRegions() const { return theJetPairRegions; }

  /**
   * Return the multi jet regions to match.
   */
  const vector<Ptr<MultiJetRegion>::ptr>& multiJetRegions() const { return theMultiJetRegions; }

  /**
   * Enumerate the ordering to apply on jets
   */
  enum Orderings {
    orderPt = 1,
    orderY = 2
  };

  /**
   * Return the ordering
   */
  int ordering() const { return theOrdering; }

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
  Ptr<MatcherBase>::ptr theUnresolvedMatcher;

  /**
   * The jet regions to match.
   */
  vector<Ptr<JetRegion>::ptr> theJetRegions;

  /**
   * The jet veto regions to check
   */
  vector<Ptr<JetRegion>::ptr> theJetVetoRegions;

  /**
   * The jet pair regions to match.
   */
  vector<Ptr<JetPairRegion>::ptr> theJetPairRegions;

  /**
   * The multi jet regions to match.
   */
  vector<Ptr<MultiJetRegion>::ptr> theMultiJetRegions;

  /**
   * Types of the ordering to apply on jets
   */
  int theOrdering;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<JetCuts> initJetCuts;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  JetCuts & operator=(const JetCuts &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of JetCuts. */
template <>
struct BaseClassTrait<JetCuts,1> {
  /** Typedef of the first base class of JetCuts. */
  typedef MultiCutBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the JetCuts class and the shared object where it is defined. */
template <>
struct ClassTraits<JetCuts>
  : public ClassTraitsBase<JetCuts> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::JetCuts"; }
  /** Return the name of the shared library be loaded to get
   *  access to the JetCuts class and every other class it uses
   *  (except the base class). */
  static string library() { return "JetCuts.so"; }
};

/** @endcond */

}

#endif /* THEPEG_JetCuts_H */
