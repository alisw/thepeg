// -*- C++ -*-
//
// InterfacedBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_InterfacedBase_H
#define ThePEG_InterfacedBase_H
// This is the declaration of the InterfacedBase class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Utilities/Named.h"
#include "ThePEG/Utilities/ClassDescription.h"
#include "ThePEG/Utilities/HoldFlag.h"
#include "InterfacedBase.xh"

namespace ThePEG {

/**
 * InterfacedBase is the base class of all Interfaced objects to be
 * handled by the BaseRepository class. InterfacedBase
 * objects can be manipulated through objects of the InterfaceBase
 * class dealing with setting parameters, switches and pointers to
 * other InterfacedBase objects.
 *
 * The InterfacedBase has a number of virtual methods to be
 * implemented by sub classes for checking the state of the object,
 * initializing the object etc.
 *
 * The InterfacedBase is derived from the PersistentBase class to
 * allow for persistent I/O, and from the Named for handling the name
 * of the object. The full name of the object is of the form
 * <code>/dir/subdir/name</code> analogous to the file name in a Unix
 * file system.
 *
 * It is possible to lock an InterfacedBase object in which case the
 * BaseRepository will not do anything that will change the state of
 * this object.
 *
 * @see BaseRepository
 * @see InterfaceBase
 */
class InterfacedBase: public PersistentBase, public Named {

  /** The BaseRepository is a close friend. */
  friend class BaseRepository;

  /** The InterfaceBase is a close friend. */
  friend class InterfaceBase;

  /** The EventGenerator is a friend. */
  friend class EventGenerator;

public:

  /**
   * Enumeration reflecting the state of an InterfacedBase object.
   */
  enum InitState {
    initializing = -1, /**< The object is currently being
			    initialized. I.e. either of update(),
			    init(), initrun() or finish() are being
			    run. */
    uninitialized = 0, /**< The object has not been initialized. */
    initialized = 1,   /**< The object has been initialized. */
    runready = 2       /**< The object is initialized and the
			    initrun() method has been called. */
  };

public:

  /**
   * The virtual (empty) destructor;
   */
  virtual ~InterfacedBase();

  /**
   * Returns the full name of this object including its path, e.g.
   * <code>/directory/subdirectory/name</code>.
   */
  string fullName() const { return Named::name(); }

  /**
   * Returns the name of this object, without the path.
   */
  string name() const {
    return Named::name().substr(Named::name().rfind('/')+1);
  }

  /**
   * Returns the path to this object including the trailing
   * '/'. <code>fullName() = path() + name()</code>.
   */
  string path() const {
    string::size_type slash = Named::name().rfind('/');
    string ret;
    if ( slash != string::npos ) ret = Named::name().substr(0,slash);
    return ret;
  }

  /**
   * Returns a comment assigned to this object.
   */
  string comment() const { return theComment; }

  /**
   * Read setup info from a standard istream \a is. May be called by
   * the Repository to initialize an object. This function first calls
   * the virtual readSetup() function to allow the sub classes the
   * part \a is to initialize themselves. What ever is left in \a is
   * after that will be assigned to the comment() of the object.
   */
  void setup(istream & is) {
    readSetup(is);
    getline(is, theComment);
  }

protected:

  /** @name Standard InterfacedBase virtual functions. */
  //@{
  /**
   * Read setup info from a standard istream \a is. May be called by
   * the Repository to initialize an object. This function is called
   * by the non virtual setup() function. A sub-class implementing it
   * should first call the base class version before parsing the \a
   * is. If the \a is is not empty after readSetup is called the
   * remaining string will be assigned to the comment() of the object.
   */
  virtual void readSetup(istream & is);

  /**
   * Check sanity of the object during the setup phase.  This function
   * is called everytime the object is changed through an interface
   * during the setup phase. Also if the setup is changed for an
   * object on which this is dependent. Note that the generator() is
   * not available when this method is called.
   *
   * This method may be called by the user interface during the setup phase
   * through the update() method after manipulating objects to check
   * the sanity of the object. When implemented by a sub class it is
   * important that the doupdate() method of the base class is called,
   * then if the sanity of this object depend on other objects, the
   * update() method of these should be called. Then if touched() is
   * true for this object or for the ones on which this depends, it is
   * an indication that some things have changed since last time
   * doupdate() was called, and the actual checking of the state of
   * this object is called for. To avoid circular loops, it is
   * important that the doupdate() method is called for the base
   * class, while the update() method is called for other objects.
   * @throws UpdateException if the setup is such that the object
   * would not work properly.
   */
  virtual void doupdate() {}

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk. Nothing should have changed since the
   * last update() call.
   * 
   * This method is called after the setup
   * phase through the init() method to indicate that the setup of a
   * run is finished. This is typpically done in a setup program
   * before this object has been saved to a run file. It must
   * therefore be made sure that the state of this object after this
   * method has been executed will not be changed if it is written to
   * a file and read in again. When implemented by a sub class it is
   * important that the doinit() method of the base class is called
   * first and then, if the initialization of this object depends on
   * other objects, that the init() method of these objects are
   * called. Only then should the class-local initialization
   * proceed. To avoid circular loops, it is important that the
   * doinit() method is called for the base class, while the init()
   * method is called for other objects.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() {}

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   *
   * This method is called just before
   * running starts through the initrun() method to indicate that the
   * actual running is to start. When implemented by a sub class it is
   * important that the doinitrun() method of the base class is called
   * first and then, if the initialization of this object depends on
   * other objects, that the initrun() method of these objects are
   * called. Only then should the class-local initialization
   * proceed. To avoid circular loops, it is important that the
   * doinitrun() method is called for the base class, while the
   * initrun() method is called for other objects.
   */
  virtual void doinitrun() {}

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   *
   * This method is called after the running
   * phase through the finish() and can eg. be used to write out
   * statistics. When implemented by a sub class it is important that
   * the dofinish() method of the base class is called while the
   * finish() methd is called for other objects.
   */
  virtual void dofinish() {}

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences() { return IVector(); }

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap &) {}
  //@}

public:

  /** @name Inlined access function. */
  //@{
  /**
   * Calls the doupdate() function with recursion prevention.
   */
  void update() {
    if ( initState ) return;
    HoldFlag<InitState> hold(initState, initializing, initialized);
    doupdate();
  }

  /**
   * Calls the doinit() function with recursion prevention.
   */
  void init() {
    if ( initState ) return;
    HoldFlag<InitState> hold(initState, initializing, initialized);
    doinit();
  }

  /**
   * Return true if this object needs to be initialized before all
   * other objects (except those for which this function also returns
   * true).  This default version always returns false, but subclasses
   * may override it to return true.
   */
  virtual bool preInitialize() const;

  /**
   * Calls the doinitrun() function with recursion prevention.
   */
  void initrun() {
    if ( initState == runready || initState == initializing ) return;
    HoldFlag<InitState> hold(initState, initializing, runready);
    doinitrun();
  }

  /**
   * Calls the dofinish() function with recursion prevention.
   */
  void finish() {
    if ( initState == uninitialized || initState == initializing ) return;
    HoldFlag<InitState> hold(initState, initializing, uninitialized);
    dofinish();
  }

  /**
   * This function should be called every time something in this
   * object has changed in a way that a sanity check with update() is
   * needed
   */
  void touch() { isTouched = true; }

  /**
   * Set the state of this object to uninitialized.
   */
  void reset() { initState = uninitialized; }

  /**
   * Calls reset() and unTouch().
   */
  void clear() {
    reset();
    untouch();
  }

  /**
   * Return the state of initialization of this object.
   */
  InitState state() const { return initState; }

  /**
   * Return true if the BaseRepository is not allowed to change the
   * state of this object.
   */
  bool locked() const { return isLocked; }

  /**
   * Return true if the state of this object has been changed since
   * the last call to update().
   */
  bool touched() const { return isTouched; }
  //@}

  /**
   * Return a full clone of this object possibly doing things to the
   * clone to make it sane.
   */
  virtual IBPtr fullclone() const { return clone(); }

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
   * Return a simple clone of this object. Should be implemented as
   * <code>return new_ptr(*this);</code> by a derived class.
   */
  virtual IBPtr clone() const = 0;

  /**
   * Protected default constructor.
   */
  InterfacedBase()
    : Named(""), isLocked(false), isTouched(true), 
      initState(uninitialized) {}

  /**
   * Protected constructor with the name given as argument.
   */
  InterfacedBase(string newName)
    : Named(newName), isLocked(false), isTouched(true),
      initState(uninitialized) {}

  /**
   * Protected copy-constructor.
   */
  InterfacedBase(const InterfacedBase & i)
    : Base(i), Named(i), isLocked(false), isTouched(true), 
      initState(uninitialized), theComment(i.theComment),
      objectDefaults(i.objectDefaults) {}

private:

  /**
   * Set a new name (full name including path).
   */
  void name(string newName) { Named::name(newName); }

  /**
   * Lock this object.
   */
  void lock() { isLocked = true; }

  /**
   * Unlock this object.
   */
  void unlock() { isLocked = false; }

  /**
   * Clear the isTouched flag.
   */
  void untouch() { isTouched = false; }

private:

  /**
   * Used by the interface to add comments.
   */
  string addComment(string);

private:

  /**
   * True if this object is not to be changed by the user interface.
   */
  bool isLocked;

  /**
   * True if this object has been chaged since the last call to
   * update().
   */
  bool isTouched;

  /**
   * Indicate if this object has been initialized or not, or if it is
   * being initialized.
   */
  InitState initState;

  /**
   * A comment assigned to this object.
   */
  string theComment;

  /**
   * A map listing object-specific defaults set for the given interfaces.
   */
  map<string,string> objectDefaults;

public:

  /**
   * Print out debugging information for this object on std::cerr. To
   * be called from within a debugger via the debug() function.
   */
  virtual void debugme() const;

private:

  /**
   * Standard Initialization object.
   */
  static AbstractClassDescription<InterfacedBase> initInterfacedBase;

  /**
   *  Private and non-existent assignment operator.
   */
  InterfacedBase & operator=(const InterfacedBase &) = delete;

protected:

  /**
   * Functor class to be used to update a range of dependent object.
   */
  struct UpdateChecker {
    /** Constructor. */
    UpdateChecker(bool & touched) : isTouched(touched) {}
    /** Constructor. */
    UpdateChecker(const UpdateChecker & uc) : isTouched(uc.isTouched) {}
    /** Call the check function for an object. */
    static void check(tIBPtr, bool &);
    /** Function call operator. */
    template <typename ptr> void operator()(const ptr & i) {
      check(i, isTouched);
    }
    /** set to false if any check() call fails. */
    bool & isTouched;
  };

  /**
   * Functor class to be used to update a range of dependent object in a map.
   */
  struct UpdateMapChecker {
    /** Constructor. */
    UpdateMapChecker(bool & touched) : isTouched(touched) {}
    /** Constructor. */
    UpdateMapChecker(const UpdateMapChecker & uc) : isTouched(uc.isTouched) {}
    /** Function call operator. */
    template <typename ref> void operator()(const ref & i) {
      UpdateChecker::check(i.second, isTouched);
    }
    /** Reference to the bool variable to be set. */
    bool & isTouched;
  };

};

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of InterfacedBase.
 */
template <>
struct BaseClassTrait<InterfacedBase,1>: public ClassTraitsType {
  /** Typedef of the base class of InterfacedBase. */
  typedef PersistentBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * InterfacedBase class.
 */
template <>
struct ClassTraits<InterfacedBase>: public ClassTraitsBase<InterfacedBase> {
  /** Return the class name. */
  static string className() { return "ThePEG::InterfacedBase"; }
};

/** @endcond */

}

#endif /* ThePEG_InterfacedBase_H */
