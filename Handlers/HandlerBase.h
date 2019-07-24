// -*- C++ -*-
//
// HandlerBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_HandlerBase_H
#define ThePEG_HandlerBase_H
// This is the declaration of the HandlerBase class.

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Repository/UseRandom.fh"
#include "ThePEG/Repository/EventGenerator.h"
#include <stdexcept>

namespace ThePEG {

template <typename T = UseRandom>
/**
 * HandlerBaseT is a dummy abstract templated class used as base class
 * to HandlerBase. HandlerBaseT inherits from the Interfaced class
 * adding some functionality such as easy acces to the RandomGenerator
 * and the StandardModel object of the controlling EventGenerator
 * object. The HandlerBaseT should only be used by the HandlerBase as
 * a base class. The fact that it is templated allows classes which in
 * turn inherits from HandlerBase to not explicitly depend on
 * EventGenerator class if the inlined accessor funtions are not
 * actually used. The only class which actually works as a template
 * argument is UseRandom, which is used to generate random numbers.
 *
 * @see Interfaced
 * @see RandomGenerator
 * @see StandardModel
 * @see EventGenerator
 * 
 */
class HandlerBaseT: public Interfaced {
public:

  /** HandlerBase is a friend. */
  friend class HandlerBase;

private:

  /** @name Standard constructors and destructors are private and can
   * only be used from the HandlerBase class. */
  //@{
  /**
   * Default constructor.
   */
  HandlerBaseT() {}

public:
  /**
   * Destructor.
   */
  virtual ~HandlerBaseT() {}
  //@}

public:

  /**
   * Return a simple flat random number in the range ]0,1[.
   */
  double rnd() const {  return T::rnd(); }

  /**
   * Return a simple flat random number in the range ]0,\a xu[.
   */
  double rnd(double xu) const { return T::rnd(xu); }

  /**
   * Return a simple flat random number in the range ]\a xl,\a xu[.
   */
  double rnd(double xl, double xu) const { return T::rnd(xl, xu); }

  /**
   * Return true with 50% probability.
   */
  bool rndbool() const { return T::rndbool(); }

  /**
   * Return a true with probability \a p.
   */
  bool rndbool(double p) const { return T::rndbool(p); }

  /**
   * Return a true with probability \a p1/(\a p1+\a p2).
   */
  bool rndbool(double p1, double p2) const { return T::rndbool(p1, p2); }

  /**
   * Return -1, 0, or 1 with relative probabilities \a p1, \a p2, \a p3.
   */
  int rndsign(double p1, double p2, double p3) const { return T::rndsign(p1, p2, p3); }

  /**
   * Return an integer \f$i\f$ with probability p\f$i\f$/(\a p0+\a p1).
   */
  int rnd2(double p0, double p1) const { return T::rnd2(p0, p1); }

  /**
   * Return an integer \f$i\f$ with probability p\f$i\f$/(\a p0+\a
   * p1+\a p2).
   */
  int rnd3(double p0, double p1, double p2) const { return T::rnd3(p0, p1, p2); }

  /**
   * Return an integer/ \f$i\f$ with probability p\f$i\f$(\a p0+\a
   * p1+\a p2+\a p3).
   */
  int rnd4(double p0, double p1, double p2, double p3) const { return T::rnd4(p0, p1, p2, p3); }

  /**
   * Return a simple flat random integrer number in the range [0,\a xu[.
   */
  long irnd(long xu = 2) const { return T::irnd(xu); }

  /**
   * Return a simple flat random integrer number in the range [\a xl,\a xu[.
   */
  long irnd(long xl, long xu) const { return T::irnd(xl, xu); }

  /**
   * Return a reference to the object containing the standard model
   * parameters for this run.
   */
  const StandardModelBase & SM() const { return *standardModel(); }

  /**
   * Return a pointer to the object containing the standard model
   * parameters for this run.
   */
  tSMPtr standardModel() const { return generator()->standardModel(); }
};

/**
 * HandlerBase is an abstract base class derived from the Interfaced
 * class via the HandlerBaseT class adding some functionality such as
 * easy acces to the RandomGenerator and the StandardModel object of
 * the controlling EventGenerator object.
 *
 * @see Interfaced
 * @see RandomGenerator
 * @see StandardModel
 * @see EventGenerator
 * 
 */
class HandlerBase: public HandlerBaseT<UseRandom> {

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

private:

  /**
   * Describe an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<HandlerBase> initHandlerBase;

  /**
   *  Private and non-existent assignment operator.
   */
  HandlerBase & operator=(const HandlerBase &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of HandlerBase.
 */
template <>
struct BaseClassTrait<HandlerBase,1>: public ClassTraitsType {
  /** Typedef of the base class of HandlerBase. Note that HandlerBaseT
   *  is not treated as a base class in this respect. */
  typedef Interfaced NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * HandlerBase class.
 */
template <>
struct ClassTraits<HandlerBase>: public ClassTraitsBase<HandlerBase> {
  /** Return the class name. */
  static string className() { return "ThePEG::HandlerBase"; }
};

/** @endcond */

}

#endif /* ThePEG_HandlerBase_H */
