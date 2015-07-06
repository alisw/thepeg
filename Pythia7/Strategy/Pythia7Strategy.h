// -*- C++ -*-
#ifndef PYTHIA7_Pythia7Strategy_H
#define PYTHIA7_Pythia7Strategy_H
// This is the declaration of the Pythia7Strategy class.

#include "ThePEG/Repository/Strategy.h"
#include "ThePEG/PDT/ParticleData.h"
// #include "Pythia7Strategy.fh"
// #include "Pythia7Strategy.xh"

namespace Pythia7 {

using namespace ThePEG;

/**
 * The Pythia7Strategy class is a sub-class of the Strategy
 * class, simply implementing the correct citation for Pythia7 in the
 * ClassDocumentation interface.
 *
 * See also \ref Pythia7StrategyInterfaces "the interfaces" defined
 * for Pythia7Strategy.
 */
class Pythia7Strategy: public Strategy {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline Pythia7Strategy();

  /**
   * Copy-constructor.
   */
  inline Pythia7Strategy(const Pythia7Strategy &);

  /**
   * Destructor.
   */
  virtual ~Pythia7Strategy();
  //@}

public:

  /**
   * Standard Init function used to initialize the interface.
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
   * Describe concrete class without persistent data.
   */
  static NoPIOClassDescription<Pythia7Strategy> initPythia7Strategy;

  /**
   *  Private and non-existent assignment operator.
   */
  Pythia7Strategy & operator=(const Pythia7Strategy &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of Pythia7::Pythia7Strategy. */
template <>
struct BaseClassTrait<Pythia7::Pythia7Strategy,1>: public ClassTraitsType {
  /** Typedef of the first base class of Pythia7::Pythia7Strategy. */
  typedef Strategy NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  Pythia7::Pythia7Strategy class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<Pythia7::Pythia7Strategy>
  : public ClassTraitsBase<Pythia7::Pythia7Strategy> {
  /** Return a platform-independent class name */
  static string className() { return "Pythia7::Pythia7Strategy"; }
  /** Return the name of the shared library be loaded to get access to
   *  the Pythia7::Pythia7Strategy class and every other class it uses
   *  (except the base class). */
  static string library() { return "Pythia7Strategy.so"; }
};

/** @endcond */

}

#include "Pythia7Strategy.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "Pythia7Strategy.tcc"
#endif

#endif /* PYTHIA7_Pythia7Strategy_H */
