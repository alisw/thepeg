// -*- C++ -*-
//
// Main.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_Main_H
#define THEPEG_Main_H
// This is the declaration of the Main class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Repository/EventGenerator.h"

namespace ThePEG {

/**
 * This is a base class for classes implementing a main steering
 * routine for running an EventGenerator, in case the standard
 * 'go()' function in the EventGenerator is not enough. An example
 * is a EventGenerator where the class deriving from Main specifies an
 * initial set of partons which should be cascaded and hadronized etc.
 *
 * The derived class will be dynamically loaded by the standard
 * <code>runThePEG</code> program and the static Init() function will
 * be run. The loaded EventGenerator will then be available in the
 * static eventGenerator() function. In addition the number of events
 * requested from the command line will be available from the N()
 * function. All other command-line arguments will be available from
 * the arguments() function.
 *
 * @see EventGenerator
 */
class Main: public Base {

public:

  /**
   * Set the global event generator.
   */
  static void eventGenerator(tEGPtr eg) { theEventGenerator = eg; }

  /**
   * Get the global event generator.
   */
  static tEGPtr eventGenerator() { return theEventGenerator; }

  /**
   * Set the number of events requested.
   */
  static void N(long n) { theN = n; }

  /**
   * Get the number of events requested.
   */
  static long N() { return theN; } 


  /**
   * Set the command-line arguments.
   */
  static void arguments(const vector<string> & args) 
  {
    theArguments = args;
  }

  /**
   * Get the command-line arguments.
   */
  static const vector<string> & arguments() { return theArguments; }

public:

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init() {}

private:

  /**
   * The global event generator.
   */
  static EGPtr theEventGenerator;

  /**
   * The number of events requested.
   */
  static long theN;

  /**
   * The command-line arguments.
   */
  static vector<string> theArguments;

private:

  /**
   * Describe an abstract base class without persistent data.
   */
  static AbstractNoPIOClassDescription<Main> initMain;

  /**
   * Private and non-existent assignment operator.
   */
  Main & operator=(const Main &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of Main. */
template <>
struct BaseClassTrait<Main,1>: public ClassTraitsType {
  /** Typedef of the first base class of Main. */
  typedef int NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  Main class. */
template <>
struct ClassTraits<Main>: public ClassTraitsBase<Main> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::Main"; }

};

/** @endcond */

}

#endif /* THEPEG_Main_H */
