// -*- C++ -*-
//
// CurrentGenerator.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_CurrentGenerator_H
#define ThePEG_CurrentGenerator_H
// This is the declaration of the CurrentGenerator class.

#include "ThePEG/Repository/EventGenerator.h"
#include "CurrentGenerator.fh"

namespace ThePEG {

/**
 * This CurrentGenerator class keeps a static stack of EventGenerators
 * which can be used anywhere by any class. When an EventGenerator is
 * initialized or run it adds itself to the stack which can be used by
 * any other object being initialized or run through the static
 * functions of the CurrentGenerator class. If someone
 * needs to use an alternative EventGenerator object a new
 * CurrentGenerator object can be constructed with a
 * pointer to the desired EventGenerator object as argument and that
 * object will the be used by the static CurrentGenerator
 * functions until the CurrentGenerator object is destructed.
 *
 * @see EventGenerator
 * 
 */
class CurrentGenerator {

public:

  /**
   * Default constructor does nothing.
   */
  CurrentGenerator() : generatorPushed(false) {}

  /**
   * Copy-constructor does nothing.
   */
  CurrentGenerator(const CurrentGenerator &)
    : generatorPushed(false) {}

  /**
   * Construct a new object specifying a new EventGenerator, \a eg, to
   * be used during this objects lifetime.
   */
  CurrentGenerator(const EGPtr & eg) : generatorPushed(false) {
    if ( eg ) {
      theGeneratorStack.push_back(eg);
      generatorPushed = true;
    }
  }

  /**
   * The destructor removing the EventGenerator specified in the
   * constructor from the stack.
   */
  ~CurrentGenerator() {
    if ( generatorPushed ) theGeneratorStack.pop_back();
  }

public:

  /**
   * Returns true if there is no currently chosen EventGenerator
   * object.
   */
  static bool isVoid() {
    return theGeneratorStack.empty() || !(theGeneratorStack.back());
  }

  /**
   * Return a reference to the currently chosen EventGenerator object.
   */
  static EventGenerator & current() {
    return *theGeneratorStack.back();
  }

  /**
   * Return a reference to the currently chosen object.
   */
  EventGenerator & operator*() const {
    return *theGeneratorStack.back();
  }

  /**
   * Return a pointer to the currently chosen object.
   */
  EventGenerator * operator->() const {
    return theGeneratorStack.back();
  }

  /**
   *  Pointer to the stack
   */
  static EventGenerator * ptr() {
    return theGeneratorStack.back();
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

  /**
   * Return a pointer to the standard model parameters used by the
   * current generator.
   */
  static tSMPtr standardModel() {
    return current().standardModel();
  }

  /**
   * Return a pointer to the strategy object containing eg. a set of
   * non-default particles to be used by the current generator.
   */
  static tStrategyPtr strategy() { return current().strategy(); }

  /**
   * Get the current standard output stream. Return a reference to the
   * stream connected to the file for general output of the current
   * generator. If no file is connected, the BaseRepository::cout()
   * will be used instead.
   */
  static ostream & out() { return current().out(); }

  /**
   * Get the current standard log stream. Return a reference to the
   * stream connected to the file for logging information of the
   * current generator. If no file is connected, the
   * BaseRepository::clog() will be used instead.
   */
  static ostream & log() { return current().log(); }

  /**
   * Get the current standard ref stream. Return a reference to the
   * stream connected to the file for references from used objects of
   * the current generator. If no file is connected, the
   * BaseRepository::cout() will be used instead.
   */
  static ostream & ref() { return current().ref(); }

  /**
   * Get object. Return a garbage collected pointer to a given object
   * in the current EventGenerator. If the object is not found, a null
   * pointer will be returned.
   */
  template <typename T>
  static typename Ptr<T>::pointer getPtr(const T & t) {
    return current().getPtr(t);
  }

  /**
   * Get object. Return a pointer to an object present in the current
   * EventGenerator given its full name. Return the null pointer if
   * non-existent.
   */
  static IBPtr getPointer(string name) {
    return current().getPointer(name);
  }

  /**
   * Get object. Return a pointer to an object of type T present in
   * the current EventGenerator given its full name. Return the null
   * pointer if non-existent.
   */
  template <typename T>
  static typename Ptr<T>::pointer getObject(string name) {
    return current().getObject<T>(name);
  }

  /**
   * Get default object. Return the default object for class T in the
   * current EventGenerator. Returns the null pointer if non-existent.
   */
  template <typename T>
  static typename Ptr<T>::pointer getDefault() {
    return current().getDefault<T>();
  }

public:

  /**
   * Class used to temporarily redirect a given ostream to the misc()
   * stream of the current EventGenerator.
   */
  class Redirect {

  public:

    /**
     * Constructor taking the stream to be redirected as input. If the
     * \a internal flag false the output will be stored in the Event
     * Generator and written to the log file in the end of the run. If
     * \internal is true the output is instead stored internally in
     * this object and is accessible through the str() function until
     * the object is destroyed.
     */
    Redirect(ostream & os, bool internal = false)
      : theStream(&os), theBuffer(os.rdbuf()) {
      if ( internal ) theStream->rdbuf(intStream.rdbuf());
      else if ( !current().useStdOut() )
	theStream->rdbuf(current().misc().rdbuf());
    }

    /**
     * The destructor which restores the original destination of the
     * stream.
     */
    ~Redirect() {
      theStream->rdbuf(theBuffer);
    }

    /**
     * If output is stored internally, acces what has been written so
     * far.
     */
    string str() const {
      return intStream.str();
    }

    /**
     * The stream which is redirected.
     */
    ostream * theStream;

    /**
     * The original buffer of the redirected stream.
     */
    std::streambuf * theBuffer;

    /**
     * An internal buffer, the content of which will be discarded when
     * the this object is destructed.
     */
    ostringstream intStream;

  };

private:

  /**
   * The stack of EventGenerators requested.
   */
  static vector<EGPtr> theGeneratorStack;

  /**
   * True if this object is responsible for pushing a EventGenerator
   * onto the stack.
   */
  bool generatorPushed;

private:

  /**
   *  Private and non-existent assignment operator.
   */
  CurrentGenerator & operator=(const CurrentGenerator &);

};

}

#endif /* ThePEG_CurrentGenerator_H */
