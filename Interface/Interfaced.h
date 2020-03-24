// -*- C++ -*-
//
// Interfaced.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Interfaced_H
#define ThePEG_Interfaced_H
// This is the declaration of the Interfaced class.

#include "ThePEG/Config/ThePEG.h"
#include "InterfacedBase.h"
#include "ThePEG/PDT/PID.h"

namespace ThePEG {

/**
 * The Interfaced class is derived from the InterfacedBase class
 * adding a couple of things particular to ThePEG, in an attempt to
 * keep the InterfacedBase class as generic as possible.
 *
 * The main addition is that the Interfaced class has a pointer to an
 * EventGenerator object. During the run-phase this points to the
 * EventGenerator controlling the run in which the Interfaced object
 * is used. Through this EventGenerator there is quick access to
 * eg. the set of ParticleData objects used, and the default
 * RandomGenerator for the run. Note that no EventGenerator object is
 * available to the Interfaced object during the setup
 * phase.
 *
 * @see InterfacedBase
 * @see EventGenerator
 * @see ParticleData
 * @see RandomGenerator
 */
class Interfaced: public InterfacedBase {

  /** Repository is a friend. */
  friend class Repository;

  /** EventGenerator is a friend. */
  friend class EventGenerator;

public:

  /**
   * Empty virtual destructor
   */
  virtual ~Interfaced();

  /**
   * Functions which are to be used during the actual event
   * generation, after the setup is complete.
   */
public:

  /**
   * A sub class can implement this function to implement some default
   * initialization for this object during the setup phase. A typical
   * example is if this object need some references to other objects
   * and if these can be easily created. In this case the objects can
   * be added to the repository in a sub-directory with the same name
   * as this object.  
   * @return false if the initialization failed.
   */
  virtual bool defaultInit();

  /** @name Functions used during the actual event generation, after
      the setup is complete. */
  //@{
  /**
   * Create a new Particle instance given a id number.
   */
  PPtr getParticle(PID) const;

  /**
   * Return a pointer to the ParticleData object corresponding to the
   * given id number.
   */
  PDPtr getParticleData(PID) const;

  /**
   * Returns true if this object has actally been used.
   */
  bool used() const { return theUseFlag; }

  /**
   * Should be called to indicate that this object has actually been
   * used.
   */
  void useMe() const { if ( !used() ) setUsed(); }

  /**
   * Return a pointer to the EventGenerator controlling the run.
   * During the setup phase this returns a null pointer.
   */
  tEGPtr generator() const { return theGenerator; }
  //@}

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
   * Standard Init function.
   */
  static void Init();

protected:

  /**
   * Register an Interfaced object with the Repository.
   */
  static void registerRepository(IBPtr);

  /**
   * Register an Interfaced object with the Repository, giving it a
   * name.
   */
  static void registerRepository(IBPtr, string newName);

  /**
   * Register the given \a object in the Repository with the given \a
   * name in a subdirectory with the same name as this object. If an
   * object with that name already exists it will be removed unless
   * there are other objects referring to it, in which case it will be
   * renamed.
   */
  void reporeg(IBPtr object, string name) const;

  /**
   * If the pointer, \a ptr, to an object is not set, create an object
   * of class \a classname and register it with the Repository with
   * the given \a objectname in a sib-directory with the same name as
   * this object.
   */
  template <typename PtrT>
  bool setDefaultReference(PtrT & ptr, string classname, string objectname) {
    if ( ptr ) return true;
    const ClassDescriptionBase * db = DescriptionList::find(classname);
    if ( !db ) return false;
    ptr = dynamic_ptr_cast<PtrT>(db->create());
    if ( !ptr ) return false;
    reporeg(ptr, objectname);
    if ( !ptr->defaultInit() ) return false;
    return true;
  }

  /**
   * Protected default constructor.
   */
  Interfaced() : theUseFlag(false) {}

  /**
   * Protected constructor taking a name as argument.
   */
  Interfaced(const string & newName) 
    : InterfacedBase(newName), theUseFlag(false) {}

  /**
   * Protected copy-constructor.
   */
  Interfaced(const Interfaced & i)
    : InterfacedBase(i), theGenerator(i.theGenerator), theUseFlag(false) {}
  
protected:

  /**
   * Protected function to reset the generator pointer, required
   * for automatic decayer generation in Herwig++ BSM models
   */
  void setGenerator(tEGPtr generator) { theGenerator=generator; }

private:

  /**
   * Used internally by 'useMe'
   */
  void setUsed() const;

  /**
   * A pointer to the EventGenerator controlling the run.
   */
  tEGPtr theGenerator;

  /**
   * Flag to tell whether this object has been used or not.
   */
  mutable bool theUseFlag;

  /**
   * Command interface function which calls defaultInit().
   */
  string doDefaultInit(string);

private:

  /**
   * Standard Initialization object.
   */
  static AbstractClassDescription<Interfaced> initInterfaced;

  /**
   *  Private and non-existent assignment operator.
   */
  Interfaced & operator=(const Interfaced &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of Interfaced.
 */
template <>
struct BaseClassTrait<Interfaced,1>: public ClassTraitsType {
  /** Typedef of the base class of Interfaced. */
  typedef InterfacedBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Interfaced class.
 */
template <>
struct ClassTraits<Interfaced>: public ClassTraitsBase<Interfaced> {
  /** Return the class name. */
  static string className() { return "ThePEG::Interfaced"; }
};

/** @endcond */

}

#endif /* ThePEG_Interfaced_H */
