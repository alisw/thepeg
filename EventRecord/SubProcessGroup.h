// -*- C++ -*-
//
// SubProcessGroup.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
// Copyright (C) 2009-2019 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_SubProcessGroup_H
#define ThePEG_SubProcessGroup_H
// This is the declaration of the SubProcessGroup class.

#include "ThePEG/EventRecord/SubProcess.h"

namespace ThePEG {

/**
 * A SubProcessGroup object represents a group of SubProcess
 * objects in dependence of a head SubProcess object.
 *
 * @see StdXCombGroup
 * @see MEGroup
 */
class SubProcessGroup: public SubProcess {

public:

  /**
   * Standard constructor.
   * @param newIncoming the two incoming partons.
   * @param newCollision the Collision to which this SubProcessGroup belongs.
   * @param newHandler the MEBase object which generated this SubProcessGroup.
   */
  SubProcessGroup(const PPair & newIncoming,
		  tCollPtr newCollision = tCollPtr(),
		  tcEventBasePtr newHandler = tcEventBasePtr());

  /**
   * Destructor.
   */
  virtual ~SubProcessGroup();

  /**
   * Return a clone of this sub process group.
   */
  virtual SubProPtr clone() const;

protected:

  /**
   * Rebind to cloned objects. When a SubProcessGroup is cloned, a shallow
   * copy is done first, then all <code>Particle</code>s etc, are
   * cloned, and finally this method is used to see to that the
   * pointers in the cloned SubProcessGroup points to the cloned
   * <code>Particle</code>s etc.
   */
  virtual void rebind(const EventTranslationMap & trans);

public:

  /**
   * Perform a LorentzTransformation of all particles in the sub
   * process.
   */
  virtual void transform(const LorentzRotation &);

  /**
   * Return the dependent SubProcess objects
   */
  const SubProcessVector& dependent() const { return theDependent; }

  /**
   * Access the dependent SubProcess objects
   */
  SubProcessVector& dependent() { return theDependent; }

  /**
   * Add a dependent SubProcess
   */
  void add(tSubProPtr sub) { dependent().push_back(sub); }

public:

  /**
   * Standard function for writing to a persistent stream.
   */
  void persistentOutput(PersistentOStream &) const;

  /**
   * Standard function for reading from a persistent stream.
   */
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function. @see Base::Init().
   */
  static void Init();

private:

  /**
   * The dependent subprocesses
   */
  SubProcessVector theDependent;

public:

  /**
   * Put to ostream
   */
  virtual void printMe(ostream&) const;

private:

  /**
   * Describe concrete class with persistent data.
   */
  static ClassDescription<SubProcessGroup> initSubProcessGroup;

  /**
   * Private default constructor must only be used by the
   * PersistentIStream class via the ClassTraits<SubProcessGroup> class .
   */
  SubProcessGroup() : SubProcess() {}

  /**
   * The ClassTraits<SubProcessGroup> class must be a friend to be able to
   * use the private default constructor.
   */
  friend struct ClassTraits<SubProcessGroup>;

  /**
   * Assignment is forbidden.
   */
  SubProcessGroup & operator=(const SubProcessGroup &) = delete;

};


/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base class of Collision. */
template <>
struct BaseClassTrait<SubProcessGroup,1>: public ClassTraitsType {
  /** Typedef of the first base class of SubProcessGroup. */
  typedef EventRecordBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SubProcessGroup class and how to create it. */
template <>
struct ClassTraits<SubProcessGroup>: public ClassTraitsBase<SubProcessGroup> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::SubProcessGroup"; }
  /** Create a SubProcessGroup object. */
  static TPtr create() { return TPtr::Create(SubProcessGroup()); }
};

/** @endcond */

}

#endif /* ThePEG_SubProcessGroup_H */
