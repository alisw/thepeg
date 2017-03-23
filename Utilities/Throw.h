// -*- C++ -*-
//
// Throw.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Throw_H
#define ThePEG_Throw_H
// This is the declaration of the Throw class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Repository/Repository.h"


namespace ThePEG {
/**
 * Helper function to make it easier to throw exceptions. The template
 * argument should be a class inheriting from Exception. In the
 * constructor, an object of this Exception class is created and
 * afterwards a message may be added using ostream-like output
 * (&lt;&lt;). If a Exception::Severity value is output with &lt;&lt;
 * the Exception object is assumed to be complete and the exception is
 * actually thrown, except if the Exception::Severity value
 * Exception::warning was specified, in which case the Exception
 * object is treated as a warning which is logged with the current
 * EventGenerator. If no current EventGenerator is present the warning
 * message is instead written to std::cerr. If no Exception::Severity
 * is specified, the Exception object is thrown when the Throw object
 * is destroyed.
 *
 * Assuming you have an Exception class called MyEx the Throw class is
 * used as follows:<br><code>Throw&lt;MyEx&gt>() &lt;&lt; "My error
 * message" &lt;&lt; Exception::eventerror;</code><br> This will throw
 * an exception and the current event will be discarded. Changing
 * <code>Exception::eventerror</code> to
 * <code>Exception::warning</code> will write out a warning, but no
 * proper exception is thrown.
 */
template <typename Ex>
struct Throw {

  /**
   * Standard constructor creating an internal Exception object.
   */
  Throw(): ex(Ex()), handled(false) {}

  /**
   * Add information to the current Exception object.
   */
  template <typename T> Throw & operator<<(const T & t) {
    ex << t;
    return *this;
  }

  /**
   * Specify the Exception::Severity of the exception. If this is
   * Exception::warning, the exception will not be thown, instead it
   * will be logged with the current * EventGenerator. If no current
   * EventGenerator is present the warning * message is instead
   * written to std::cerr. All other seveities will cause the
   * exception to be thrown immediately.
   */
  void operator<<(Exception::Severity sev) {
    handled = true;
    ex << sev;
    if ( sev != Exception::warning && sev != Exception::info  ) {
      throw ex;
    } else {
      if ( CurrentGenerator::isVoid() ) {
	Repository::clog() << ex.message() << endl;
	ex.handle();
      } else {
	CurrentGenerator::current().logWarning(ex);
      }
    }
  }

  /**
   * The destructor will throw the exception if it has not been handled.
   */
  ~Throw() throw (Ex) {
    if ( !handled ) throw ex;
  }

  /**
   * The ExceptionObject to be thrown.
   */
  Ex ex;

  /**
   * If true, the exception has been handled and should not be thrown
   * in the destructor.
   */
  bool handled;
};


}

#endif /* ThePEG_Throw_H */
