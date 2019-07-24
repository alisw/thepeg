// -*- C++ -*-
//
// Direction.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Direction_H
#define ThePEG_Direction_H
// This is the declaration of the Direction class.

#include "ThePEG/Config/ThePEG.h"
#include "Direction.xh"

namespace ThePEG {

template <int I>
/**
 * A <code>Direction</code> object can be used to specify that some
 * following operations should be assumed to be performed with the
 * z-direction of the momenta reversed. As an example, if
 * <code>Direction<0>::pos()</code> is true, the method
 * <code>Lorentz5Momentum::dirPlus()</code> will return the positive,
 * light-cone component, and <code>Lorentz5Momentum::dirMinus()</code>
 * the negative, while if <code>Direction<0>::pos()</code> is false
 * the behavior of the methods are reversed.
 *
 * <code>Direction</code> is templated with an integer template argument
 * (default = 0), and only one object per class can be instatiated at
 * the time. Attempts to instatiate a second object of a
 * <code>Direction</code> class will result in an exception being
 * thrown. To have several different directions classes with different
 * template arguments must be instantiated. <code>Direction<0></code> is
 * reserved for <code>Lorentz5Momentum</code>. Attempts to use the
 * static methods of a <code>Direction</code> class when no object has
 * been instatiated will result in an exception being thrown.
 *
 * @see Lorentz5Momentum
 */
class Direction {

public:

  /** The enum defining the directions. */
  enum Dir { Neg = -1, /**< Reversed direction. */
	     Negative = -1, /**< Reversed direction. */
	     Undefined = 0, /**< No direction has been defined. */
	     Pos = 1, /**< Standard (positive) direction. */
	     Positive = 1 /**< Standard (positive) direction. */
  };

public:

  /**
   * Create an object with a given direction.
   */
  Direction(Dir newDirection)
    
  {
    if ( theDirection != Undefined ) throw MultipleDirectionException(I);
    if ( newDirection == Positive ) theDirection = Positive;
    else if ( newDirection == Negative ) theDirection = Negative;
    else throw UndefinedDirectionException(I);
  }

  /**
   * Create an object with a positive direction if rnd > 0.5,
   * otherwise set the negative direction.
   */
  Direction(double rnd)
  {
    if ( theDirection != Undefined ) throw MultipleDirectionException(I);
    theDirection = rnd > 0 ? Positive : Negative;
  }

  /**
   * Create an object with a positive direction if p is true,
   * otherwise set the negative direction.
   */
  Direction(bool p)
  {
    if ( theDirection != Undefined ) throw MultipleDirectionException(I);
    theDirection = p ? Positive : Negative;
  }

  /**
   * Destructure makeing the static variable undefined.
   */
  ~Direction() { theDirection = Undefined; }

public:

  /**
   * Set the direction.
   */
  static void set(Dir newDirection) {
    if ( newDirection == Positive ) theDirection = Positive;
    else if ( newDirection == Negative ) theDirection = Negative;
    else throw UndefinedDirectionException(I);
  }

  /**
   * Reverse the direction.
   */
  static void reverse() {
    theDirection = pos() ? Negative : Positive;
  }

  /**
   * Return true if the direction is positive.
   */
  static bool pos() {
    return dir() == Positive;
  }

  /**
   * Return true if the direction is negative (reversed).
   */
  static bool neg() {
    return dir() == Negative;
  }

  /**
   * Return the direction.
   */
  static Dir dir() {
    if ( theDirection == Undefined ) throw UndefinedDirectionException(I);
    return theDirection;
  }

private:

  /**
   * The direction.
   */
  static Dir theDirection;

private:

  /**
   * Default ctors and assignment is private and not implemented.
   */
  Direction();
  /**
   * Default ctors and assignment is private and not implemented.
   */
  Direction(const Direction &);
  /**
   * Default ctors and assignment is private and not implemented.
   */
  Direction & operator=(const Direction &) = delete;

};

template<int I>
typename Direction<I>::Dir Direction<I>::theDirection = Direction<I>::Undefined;

}

#endif /* ThePEG_Direction_H */
