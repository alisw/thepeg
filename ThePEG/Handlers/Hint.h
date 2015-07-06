// -*- C++ -*-
//
// Hint.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Hint_H
#define ThePEG_Hint_H
// This is the declaration of the Hint class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Utilities/ClassDescription.h"

namespace ThePEG {

/**
 * Hint is a base class to be used to pass information
 * between <code>StepHandler</code> s, which cannot be convayed
 * through the Event record. The base class contains a vector of of
 * tagged particles. A StepHandler is always given a hint, and is only
 * allowed to treat Particles from the current Step which are listed
 * in the vector of tagged particles in the hint (if this vector is
 * empty the StepHandler may treat all particles in the Step.
 *
 * A Hint may have the stop flag set. In this case
 * the StepHandler to which the hint is assigned is
 * not called, and the event generation is stopped.
 *
 * A Hint may be given a scale, but what a StepHandler does with this
 * and other pieces of information possibly supplied by subclasses of
 * Hint, is not defined.
 *
 * There is a special Hint which is kept as the static member called
 * Hint::theDefaultHint. Although any default constructed Hint object
 * would work as a default hint, only pointers to this static object
 * should be used where a default hint is needed.
 *
 *
 * @see StepHandler
 * @see EventHandler
 * @see Particle
 * @see Event
 * @see Step
 * 
 */
class Hint: public Base {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  Hint() : theScale(Energy2()), theStopFlag(false) {}
  //@}

public:

  /**
   * Return true if there are tagged particles in the hint.
   */
  bool tagged() const { return !theTagged.empty(); }

  /**
   * Return a list of pointers to particles to be handled. A handler
   * is not allowed to touch other particles in the event record. If a
   * particle which has been flagged by the hint is no longer present
   * in the current Step, a null pointer is inserted in its place.
   */
  tPVector tagged(const Step & s) const;

  /**
   * Add a range of particles to the list of tagged particles.
   */
  template <typename InputIterator>
  void tag(InputIterator first, InputIterator last) 
  { 
    theTagged.insert(theTagged.end(), first, last); 
  }

  /**
   * Add a particle to the list of tagged particles.
   */
  void tag(tPPtr p) { if (p) theTagged.push_back(p); }

  /**
   * Set the stop hint.
   */
  void stop(bool newStopFlag) 
  {
    theStopFlag = newStopFlag;
    if ( theStopFlag ) theTagged.clear(); 
  }

  /**
   * Get the stop hint.
   */
  bool stop() const { return theStopFlag; }

  /**
   * Set the scale.
   */
  void scale(const Scale & newScale) { theScale = newScale; }
  /**
   * Get the scale.
   */
  const Scale & scale() const { return theScale; }

  /**
   * Return a pointer to the default hint.
   */
  static tHintPtr Default() { return tHintPtr(&theDefaultHint); }

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
   * Standard Init function used to initialize the interface.
   */
  static void Init();

private:

  /**
   * The vector of tagged particles.
   */
  tPVector theTagged;

  /**
   * The scale.
   */
  Scale theScale;

  /**
   * The stop hint.
   */
  bool theStopFlag;

  /**
   * The default hint.
   */
  static Hint theDefaultHint;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<Hint> initHint;

  /**
   * Assignment is private and non-existing.
   */
  Hint & operator=(const Hint & h);

};


/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of Hint.
 */
template <>
struct BaseClassTrait<Hint,1>: public ClassTraitsType {
  /** Typedef of the base class of Hint. */
  typedef Base NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Hint class.
 */
template <>
struct ClassTraits<Hint>:
    public ClassTraitsBase<Hint> {
  /** Return the class name. */
  static string className() { return "ThePEG::Hint"; }
};

/** @endcond */

}

#endif /* ThePEG_Hint_H */
