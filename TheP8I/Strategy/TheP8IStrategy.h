// -*- C++ -*-
#ifndef THEP8I_TheP8IStrategy_H
#define THEP8I_TheP8IStrategy_H
// This is the declaration of the TheP8IStrategy class.

#include "ThePEG/Repository/Strategy.h"
#include "ThePEG/PDT/ParticleData.h"
// #include "TheP8IStrategy.fh"
// #include "TheP8IStrategy.xh"

namespace TheP8I {

using namespace ThePEG;

/**
 * The TheP8IStrategy class is a sub-class of the Strategy
 * class, simply implementing the correct citation for TheP8I in the
 * ClassDocumentation interface.
 *
 * See also \ref TheP8IStrategyInterfaces "the interfaces" defined
 * for TheP8IStrategy.
 */
class TheP8IStrategy: public Strategy {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline TheP8IStrategy() {}

  /**
   * Destructor.
   */
  virtual ~TheP8IStrategy();
  //@}

public:

  /**
   * Standard Init function used to initialize the interface.
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
   * Describe concrete class without persistent data.
   */
  static NoPIOClassDescription<TheP8IStrategy> initTheP8IStrategy;

  /**
   *  Private and non-existent assignment operator.
   */
  TheP8IStrategy & operator=(const TheP8IStrategy &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of TheP8I::TheP8IStrategy. */
template <>
struct BaseClassTrait<TheP8I::TheP8IStrategy,1>: public ClassTraitsType {
  /** Typedef of the first base class of TheP8I::TheP8IStrategy. */
  typedef Strategy NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  TheP8I::TheP8IStrategy class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<TheP8I::TheP8IStrategy>
  : public ClassTraitsBase<TheP8I::TheP8IStrategy> {
  /** Return a platform-independent class name */
  static string className() { return "TheP8I::TheP8IStrategy"; }
  /** Return the name of the shared library be loaded to get access to
   *  the TheP8I::TheP8IStrategy class and every other class it uses
   *  (except the base class). */
  static string library() { return "libTheP8I.so"; }
};

/** @endcond */

}

#endif /* THEP8I_TheP8IStrategy_H */
