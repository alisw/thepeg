// -*- C++ -*-
//
// Exception.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Exception_H
#define ThePEG_Exception_H
// This is the declaration of the Exception class.

#include <exception>
#include "ThePEG/Config/ThePEG.h"
#include "Exception.fh"
#include <string>
#include <iosfwd>

extern "C" {
  void breakThePEG();
}

namespace ThePEG {

/**
 * <code>Exception</code> is the base class for all exceptions to be
 * used in ThePEG. It is derived from <code>std::exception</code> and
 * adds information about the severity of the exception to indicate to
 * the Repository and EventGenrator how to act on it.
 *
 * To throw an exception one should inherit from
 * <code>Exception</code> and add information in the constructor of
 * the base class. Alternatively one can use the <code>operator<<
 * </code> operator on a default constructed <code>Exception</code> to
 * add information as for a standard <code>ostream</code> object, in
 * which case one should always end with adding an enum of the type
 * <code>Exception::Severity</code> to indicate the severity of the
 * exception e.g.<br> <code>Exception() << "Something went wrong." <<
 * Exception::eventerror</code>.
 *
 * @see Repository
 * @see EventGenrator
 */
class Exception: public exception {

public:

  /**
   * The levels of severity.
   */
  enum Severity {
    unknown,      /**< Unknown severity */
    info,         /**< Not severe (but the user should be
		   *   informed). */
    warning,      /**< Possibly severe, (the user should be
		   *   warned). */
    setuperror,   /**< Command failed during setup phase,
		   *   execution is continued. */
    eventerror,   /**< Possibly severe, (the event being
		   *   generated should be discarded). */
    runerror,     /**< Severe error, (the run should be
		   *   terminated). */
    maybeabort,   /**< Severe error, (the run should be
		   *   terminated, possibly dumping core). */
    abortnow      /**< Severe error, (the run is aborted
		   *   immediately, before the exception is
		   *   thrown). */
  };

public:

  /**
   * Standard constructor.
   * @param str an error message.
   * @param sev the severity.
   */
  Exception(const string & str, Severity sev);

  /**
   * Default constructor.
   */
  Exception() : handled(false), theSeverity(unknown) { breakThePEG(); }

  /**
   * The copy constructor.
   */
  Exception(const Exception & ex) 
    : std::exception(ex), theMessage(ex.message()), 
      handled(ex.handled), theSeverity(ex.severity()) 
  {
    ex.handle();
  }

  /**
   * The destructor
   */
  virtual ~Exception() noexcept;

public:

  /**
   * Assignment.
   */
  const Exception & operator=(const Exception & ex) {
    handled = ex.handled;
    theMessage << ex.message();
    theSeverity = ex.severity();
    ex.handle();
    return *this;
  }

  /**
   * Comparison
   */
  bool operator==(const Exception & ex) const {
    return ( message() == ex.message() && severity() == ex.severity() );
  }

  /**
   * Compare severity. If equal compare error message
   * lexicographically.
   */
  bool operator<(const Exception & ex) const {
    return ( severity() == ex.severity() ? 
	     ( message() < ex.message() ) :
	     ( severity() < ex.severity() ) );
  }

public:

  /**
   * Return the error message.
   */
  virtual const char* what() const noexcept {
    static string str;
    str = message();
    return str.c_str();
  }

  /**
   * Return the error message.
   */
  string message() const {
    string mess = theMessage.str();
    return mess.empty() ? string("Error message not provided.") : mess;
  }

  /**
   * Write the error message to a stream.
   */
  void writeMessage(ostream & os = *errstream) const;

  /**
   * Return the severity.
   */
  Severity severity() const { return theSeverity; }

  /**
   * Indicate that this exception has been taken care of.
   */
  void handle() const { handled = true; }

  /**
   * Add info to the exception message.
   */
  template <typename T>
  Exception & operator<<(const T & t) {
    theMessage << t;
    return *this;
  }

  /**
   * Set the severity for the exception.
   */
  Exception & operator<<(Severity sev) {
    severity(sev);
    return *this;
  }

protected:

  /**
   * set the severity.
   */
  void severity(Severity);

  /**
   * Stream to write the error message to.
   */
  mutable ostringstream theMessage;

private:

  /**
   * True if this exception has been taken care of.
   */
  mutable bool handled;

  /**
   * The severity.
   */
  Severity theSeverity;

  /**
   * The default stream to write the error message if unhandled.
   */
  static ostream * errstream;

public:

  /**
   * If this flag is set, all abortnow and maybeabort severities will
   * be treated as runerror.
   */
  static bool noabort;

};

}

#endif /* ThePEG_Exception_H */
