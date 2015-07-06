// -*- C++ -*-
#ifndef PYTHIA7_StringCollapser_H
#define PYTHIA7_StringCollapser_H
// This is the declaration of the StringCollapser class.

#include "Pythia7/Config/Pythia7.h"
#include "ThePEG/Handlers/ClusterCollapser.h"
// #include "StringCollapser.fh"
// #include "StringCollapser.xh"

namespace Pythia7 {

/**
 * StringCollapser is the class used by the LundFragHandler class to
 * collapse strings which are deemed too small to fragment into one or
 * two hadrons. Currently this class does not introduce any
 * functionality on top the ClusterCollapser base class.
 *
 * @see \ref StringCollapserInterfaces "The interfaces"
 * defined for StringCollapser.
 */
class StringCollapser: public ClusterCollapser {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline StringCollapser();

  /**
   * Copy-constructor.
   */
  inline StringCollapser(const StringCollapser &);

  /**
   * Destructor.
   */
  virtual ~StringCollapser();
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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:


protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<StringCollapser> initStringCollapser;

  /**
   *  Private and non-existent assignment operator.
   */
  StringCollapser & operator=(const StringCollapser &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * Pythia7::StringCollapser.
 */
template <>
struct BaseClassTrait<Pythia7::StringCollapser,1>: public ClassTraitsType {
  /** Typedef of the base class of Pythia7::StringCollapser. */
  typedef ClusterCollapser NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Pythia7::StringCollapser class and the shared object where it
 * is defined.
 */
template <>
struct ClassTraits<Pythia7::StringCollapser>
  : public ClassTraitsBase<Pythia7::StringCollapser> {
  /** Return the class name.  */
  static string className() { return "Pythia7::StringCollapser"; }
  /**
   * Return the name of the shared library to be loaded to get access
   * to the Pythia7::StringCollapser class and every other class
   * it uses (except the base class).
   */
  static string library() { return "libP7String.so"; }
};

/** @endcond */

}

#include "StringCollapser.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "StringCollapser.tcc"
#endif

#endif /* PYTHIA7_StringCollapser_H */
