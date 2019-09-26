// -*- C++ -*-
//
// MultiEventGenerator.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_MultiEventGenerator_H
#define ThePEG_MultiEventGenerator_H
// This is the declaration of the MultiEventGenerator class.

#include "EventGenerator.h"

namespace ThePEG {

/**
 * The MultiEventGenerator class is derived from the
 * EventGenerator class and is capable of making several runs with
 * a pre-defined set of parameter and switch values.
 *
 * With the Command<MultiEventGenerator> interface AddInterface set of
 * parameters for an included object can be specified as eg.
 * <code>object:interface arg1, arg2, arg3 ...</code>. The event
 * generator will then be run once with the specified objects
 * interface set to <code>arg1</code>, then once with
 * <code>arg2</code> etc. If several AddInterface commands are given
 * the event generator will be run once for each possible combination
 * arguments to object interfaces.
 *
 * @see EventGenerator
 * 
 */
class MultiEventGenerator: public EventGenerator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  MultiEventGenerator(): firstSubrun(0), lastSubrun(0) {}

  /**
   * Destructor.
   */
  virtual ~MultiEventGenerator();
  //@}

  /**
   * Append a tag to the run name. Derived classes may put special
   * meaning to the tags. In this case a tag on the form #n will
   * select only subrun number n, while #n-m will select subruns n
   * through m.
   */
  virtual void addTag(string tag);

protected:

  /** @name Public virtual functions required by the base class. */
  //@{
  /**
   * Run this EventGenerator session. Is called from
   * EventGenerator::go(long,long,bool).
   */
  virtual void doGo(long next, long maxevent, bool tics);
  //@}

  /** @name Functions used by the Command<MultiEventGenerator>
      interfaces to set up the different parameters of the runs. */
  //@{
  /**
   * Used to add an interface of an object which should be used with a
   * set of different values. The argument should be given as
   * "object:interface arg1, arg2, ..."
   */
  string addInterface(string);

  /**
   * Used to add an interface of an object a random value each
   * run. The argument should be given as "object:interface N min max
   * mean width"
   */
  string addRndInterface(string);

  /**
   * Helper function for addInterface(string) and addRndInterface(string).
   */
  string addInterface(string, bool rnd);


  /**
   * Used to remove an interface of an object which should be used
   * with a set of different values. The argument should be given as
   * "object:interface arg1, arg2, ..."
   */
  string removeInterface(string);
  //@}

  /**
   * return a header for this sub-run.
   */
  string heading(long, const vector<const InterfaceBase *> &, string) const;

  /**
   * A separate random number generator to be used for generating
   * random parameter values (to ensure reproducable sequences).
   */
  RandomGenerator & randomArg() const {
    return theSeparateRandom? *theSeparateRandom: random();
  }

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
   * Standard Init function used to initialize the interfaces.
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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}

private:

  /** A vector of character strings. */
  typedef vector<string> StringVector;

  /**
   * The objects for which there are different interface settings.
   */
  IVector theObjects;

  /**
   * The interfaces to be modified for the corresponding objects in
   * theObjects.
   */
  StringVector theInterfaces;

  /**
   * If the there are positional arguments to theInterfaces these are
   * specified here.
   */
  StringVector thePosArgs;

  /**
   * The arguments to be used for each of theInterfaces.
   */
  vector<StringVector> theValues;

  /**
   * If non zero, the first subrun to be considered.
   */
  int firstSubrun;

  /**
   * If non zero, the last subrun to be considered.
   */
  int lastSubrun;

  /**
   * A separate random number generator to be used for generating
   * random parameter values (to ensure reproducable sequences).
   */
  RanGenPtr theSeparateRandom;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<MultiEventGenerator> initMultiEventGenerator;

  /**
   *  Private and non-existent assignment operator.
   */
  MultiEventGenerator & operator=(const MultiEventGenerator &) = delete;

};


/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of MultiEventGenerator. */
template <>
struct BaseClassTrait<MultiEventGenerator,1>: public ClassTraitsType {
  /** Typedef of the first base class of MultiEventGenerator. */
  typedef EventGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  MultiEventGenerator class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<MultiEventGenerator>:
    public ClassTraitsBase<MultiEventGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::MultiEventGenerator"; }
  /** Return the name of the shared library be loaded to get access to
   *  the MultiEventGenerator class and every other class it uses
   *  (except the base class). */
  static string library() { return "MultiEventGenerator.so"; }

};

/** @endcond */

}

#endif /* ThePEG_MultiEventGenerator_H */
