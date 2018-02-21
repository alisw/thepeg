// -*- C++ -*-
//
// UseRandom.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_UseRandom_H
#define ThePEG_UseRandom_H
// This is the declaration of the UseRandom class.

#include "ThePEG/Repository/RandomGenerator.h"

namespace ThePEG {

/**
 * This UseRandom class keeps a static stack of RandomGenerator
 * objects which can be used anywhere by any class. When an
 * EventGenerator is initialized or run it adds a RandomGenerator
 * object to the stack which can be used by any other object being
 * initialized or run through the static functions of the UseRandom
 * class. If someone needs to use an alternative RandomGenerator
 * object a new UseRandom object can be constructed with a pointer to
 * the desired RandomGenerator object as argument and that object will
 * the be used by the static UseRandom functions until the UseRandom
 * object is destructed.
 *
 * @see RandomGenerator
 * @see EventGenerator
 * 
 */
class UseRandom {

public:

  /**
   * Default constructor does nothing.
   */
  UseRandom() : randomPushed(false) {}

  /**
   * Copy-constructor does nothing.
   */
  UseRandom(const UseRandom &) : randomPushed(false) {}

  /**
   * Construct a new object specifying a new RandomGenerator, \a r, to
   * be used during this objects lifetime
   */
  UseRandom(const RanGenPtr & r) : randomPushed(false) {
    if ( r ) {
      theRandomStack.push_back(r);
      randomPushed = true;
    }
  }

  /**
   * The destructor removing the RandomGenerator specified in the
   * constructor from the stack.
   */
  ~UseRandom() { if ( randomPushed ) theRandomStack.pop_back(); }

public:

  /**
   * Return a reference to the currently chosen RandomGenerator object.
   */
  static RandomGenerator & current() { return *theRandomStack.back(); }

  /**
   * Return a pointer to the currently chosen RandomGenerator object.
   */
//  static RandomEngine * currentEngine() {
//    return &(current().randomGenerator());
//  }

  /**
   * Return a simple flat random number (from the current
   * RandomGenerator object) in the range ]0,1[.
   */
  static double rnd() { return current().rnd(); }

  /**
   * Return \a n simple flat random number (from the current
   * RandomGenerator object) in the range ]0,1[.
   */
  static RandomGenerator::RndVector rndvec(int n) {
    return current().rndvec(n);
  }

  /**
   * Return a simple flat random number (from the current
   * RandomGenerator object) in the range ]0,\a xu[.
   */
  template <typename Unit>
  static Unit rnd(Unit xu) { return current().rnd(xu); }

  /**
   * Return a simple flat random number (from the current
   * RandomGenerator object) in the range ]\a xl,\a xu[.
   */
  template <typename Unit>
  static Unit rnd(Unit xl, Unit xu) { 
    return current().rnd(xl, xu); 
  }
  
  /**
   * Return a true with probability \a p (default 0.5).
   */
  static bool rndbool(double p = 0.5) {
    return current().rndbool(p);
  }

  /**
   * Return a true with probability \a p (default 0.5). Uses push_back
   * to reuse random number.
   */
  static bool prndbool(double p = 0.5) {
    return current().rndbool(p);
  }

  /**
   * Return a true with probability \a p1/(\a p1+\a p2).
   */
  static bool rndbool(double p1, double p2) {
    return current().rndbool(p1, p2);
  }

  /**
   * Return a true with probability \a p1/(\a p1+\a p2). Uses
   * push_back to reuse random number.
   */
  static bool prndbool(double p1, double p2) {
    return current().rndbool(p1, p2);
  }

  /**
   * Return -1, 0, or 1 with relative probabilities \a p1, \a p2, \a
   * p3.
   */
  static int rndsign(double p1, double p2, double p3) {
    return current().rndsign(p1, p2, p3);
  }

  /**
   * Return -1, 0, or 1 with relative probabilities \a p1, \a p2, \a
   * p3. Uses push_back to reuse random number.
   */
  static int prndsign(double p1, double p2, double p3) {
    return current().rndsign(p1, p2, p3);
  }

  /**
   * Return an integer \f$i\f$ with probability p\f$i\f$/(\a p0+\a
   * p1).
   */
  static int rnd2(double p0, double p1) {
    return current().rnd2(p0, p1);
  }

  /**
   * Return an integer \f$i\f$ with probability p\f$i\f$/(\a p0+\a
   * p1+\a p2).
   */
  static int rnd3(double p0, double p1, double p2) {
    return current().rnd3(p0, p1, p2);
  }

  /**
   * Return an integer/ \f$i\f$ with probability p\f$i\f$(\a p0+\a
   * p1+\a p2+\a p3).
   */
  static int rnd4(double p0, double p1, double p2, double p3) {
    return current().rnd4(p0, p1, p2, p3);
  }

  /**
   * Return a simple flat random integrer number in the range [0,\a xu[.
   */
  static long irnd(long xu = 2) { return long(rnd() * xu); }

  /**
   * Return a simple flat random integrer number in the range [\a xl,\a xu[.
   */
  static long irnd(long xl, long xu) { return xl + irnd(xu-xl); }
  
  /**
   * Return a number between zero and infinity, distributed according
   * to \f$e^-x\f$.
   */
  static double rndExp() { return current().rndExp(); }

  /**
   * Return a number between zero and infinity, distributed according
   * to \f$e^-{x/\mu}\f$ where \f$\mu\f$ is the \a mean value.
   */
  template <typename Unit>
  static Unit rndExp(Unit mean) { return current().rndExp(mean); }

  /**
   * Return a number distributed according to a Gaussian distribution
   * with zero mean and unit variance.
   */
  static double rndGauss() { return current().rndGauss(); }

  /**
   * Return a number distributed according to a Gaussian distribution
   * with a given standard deviation, \a sigma, and a given \a mean.
   */
  template <typename Unit>
  static Unit rndGauss(Unit sigma, Unit mean = Unit()) {
    return current().rndGauss(sigma, mean);
  }

  /**
   * Return a positive number distributed according to a
   * non-relativistic Breit-Wigner with a given width, \a gamma, and a
   * given \a mean.
   */
  template <typename Unit>
  static Unit rndBW(Unit mean, Unit gamma) {
    return current().rndBW(mean, gamma);
  }

  /**
   * Return a positive number distributed according to a
   * non-relativistic Breit-Wigner with a given width, \a gamma, and a
   * given \a mean. The distribution is cut-off so that the number is
   * between \a mean - \a cut and \a mean + \a cut
   */
  template <typename Unit>
  static Unit rndBW(Unit mean, Unit gamma, Unit cut) {
    return current().rndBW(mean, gamma, cut);
  }

  /**
   * Return a positive number distributed according to a relativistic
   * Breit-Wigner with a given width, \a gamma, and a given \a mean.
   */
  template <typename Unit>
  static Unit rndRelBW(Unit mean, Unit gamma) {
    return current().rndRelBW(mean, gamma);
  }

  /**
   * Return a positive number distributed according to a relativistic
   * Breit-Wigner with a given width, \a gamma, and a given \a
   * mean. The distribution is cut-off so that the number is between
   * \a mean - \a cut and \a mean + \a cut
   */
  template <typename Unit>
  static Unit rndRelBW(Unit mean, Unit gamma, Unit cut) {
    return current().rndRelBW(mean, gamma, cut);
  }

  /**
   * Return a non-negative number generated according to a Poissonian
   * distribution with a given \a mean.
   */
  static long rndPoisson(double mean) {
    return current().rndPoisson(mean);
  }

private:

  /**
   * The stack of RandomGenerators requested.
   */
  static vector<RanGenPtr> theRandomStack;

  /**
   * True if this object is responsible for pushing a RandomGenerator
   * onto the stack.
   */
  bool randomPushed;

private:

  /**
   *  Private and non-existent assignment operator.
   */
  UseRandom & operator=(const UseRandom &);

};

}

#endif /* ThePEG_UseRandom_H */
