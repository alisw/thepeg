// -*- C++ -*-
//
// Current.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Current_H
#define ThePEG_Current_H
// This is the declaration of the Current class.

namespace ThePEG {

/**
 * The Current class keeps a static stack of objects of the templated
 * class, which can be used anywhere by any class. When an object is
 * active it adds itself to the stack which can be used by any other
 * object through the static functions of the Current class. If
 * someone needs to use an alternative object a new Current object can
 * be constructed with a pointer to the desired object
 * as argument and that object will the be used by the static Current
 * functions until the Current object is destructed.
 *
 * Default-contructed objects of the Current class can be used as a
 * pointer to the currently chosen object on the stack.
 *
 * The typical use case for this class is a handler class which uses a
 * number of objects which do not have a reference back to the
 * handler, but still need to acces some member functions. In a member
 * function the handler class will construct a Current object:
 * <code>Current&lt;Handler&gt; current(this);</code> in any following
 * function called in this member function, any object can then access
 * the handlers methods as
 * <code>Current&lt;Handler&gt;()-&gt;memfun();</code>.
 *
 */
template <typename T>
class Current {

public:

  /**
   * Default constructor does nothing.
   */
  Current() : pushed(false) {}

  /**
   * Copy-constructor does nothing.
   */
  Current(const Current<T> &)
    : pushed(false) {}

  /**
   * Construct a new object specifying a new object, \a o, to be used
   * during this objects lifetime. The object must not be deleted
   * until the Current object us destroyed.
   */
  Current(T * t) : pushed(false) {
    if ( t ) {
      theStack.push_back(t);
      pushed = true;
    }
  }

  /**
   * The destructor removing the object specified in the constructor
   * from the stack.
   */
  ~Current() {
    if ( pushed ) theStack.pop_back();
  }

public:

  /**
   * Returns true if there is no currently chosen object.
   */
  static bool isVoid() {
    return theStack.empty() || !(theStack.back());
  }

  /**
   * Return a reference to the currently chosen object.
   */
  static T & current() {
    return *theStack.back();
  }

  /**
   * Return a reference to the currently chosen object.
   */
  T & operator*() const {
    return *theStack.back();
  }

  /**
   * Return a pointer to the currently chosen object.
   */
  T * operator->() const {
    return theStack.back();
  }

  /**
   *  Pointer to the stack
   */
  static T * ptr() {
    return theStack.back();
  }

  /**
   * Test for existance
   */
  operator bool() const {
    return ptr();
  }

  /**
   * Test for existance
   */
  bool operator!() const {
    return !ptr();
  }

private:

  /**
   * The stack of objects requested.
   */
  static vector<T *> theStack;

  /**
   * True if this object is responsible for pushing an object
   * onto the stack.
   */
  bool pushed;

private:

  /**
   *  Private and non-existent assignment operator.
   */
  Current<T> & operator=(const Current<T> &);

};

template <typename T>
std::vector<T *> ThePEG::Current<T>::theStack;

}

#endif /* ThePEG_Current_H */
