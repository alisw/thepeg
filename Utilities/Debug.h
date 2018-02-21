// -*- C++ -*-
//
// Debug.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Debug_H
#define ThePEG_Debug_H
// This is the declaration of the Debug class.

#include <vector>

namespace ThePEG {

/**
 * The <code>Debug</code> class defines a number of debug levels and a
 * static public member giving the current debug level in a run.
 */
class Debug {

public:

  /**
   * The different debug levels.
   */
  enum Levels {
    noDebug = 0, /**< No debugging. */
    printSomeEvents = 1, /**< Lowest debug level. Some events are printed out. */
    printEveryEvent = 5, /**< Higher debug level. All events are printed out. */
    full = 9 /**< Highest possible debug level */
  };

  /**
   * The current level.
   */
  static int level;

  /**
   * If true, the debug level has been set from the outside from the
   * calling program. This would then override any debug settings in
   * the event generator.
   */
  static bool isset;

  /**
   * A vector of switches indicating whether a given <i>debug item</i>
   * is switched on or not. The index of a debug item has no special
   * meaning. An implementor may assume that a given number
   * corresponds to a certain request for debug output, but cannot be
   * sure that someone else uses that number for some other purpose.
   */
  static std::vector<bool> debugItems;

  /**
   * Switch on or off a given debug item. If no such item exists, one
   * will be created.
   */
  static void debugItem(int item, bool on);

  /**
   * Set the debug level.
   */
  static void setDebug(int ilev);

  /**
   * Check if a given item should be debugged. If no such item is
   * present false is returned.
   */
  static bool debugItem(int item)
  {
    if ( level == noDebug ) return false;
    if ( level == full ) return true;
    return ( item < 0 || std::size_t(item) >= debugItems.size() )? false:
      debugItems[item];
  }

  /** @name Functions for manipulating floating point (FPU) signals
   * (requires the existence of <code>fpu_controll.h</code> on the
   * platform). */
  //@{
  /**
   * Switch on signalling all FPU exceptions.
   */
  static void unmaskFpuErrors();
  /**
   * Switch on signalling overflow FPU exceptions.
   */
  static void unmaskFpuOverflow();
  /**
   * Switch on signalling underflow FPU exceptions.
   */
  static void unmaskFpuUnderflow();
  /**
   * Switch on signalling division by zero FPU exceptions.
   */
  static void unmaskFpuDivZero();
  /**
   * Switch on signalling denormalized argument FPU exceptions.
   */
  static void unmaskFpuDenorm();
  /**
   * Switch on signalling invalid operation FPU exceptions.
   */
  static void unmaskFpuInvalid();
  /**
   * Switch off signalling all FPU exceptions.
   */
  static void maskFpuErrors();
  /**
   * Switch off signalling overflow FPU exceptions.
   */
  static void maskFpuOverflow();
  /**
   * Switch off signalling underflow FPU exceptions.
   */
  static void maskFpuUnderflow();
  /**
   * Switch off signalling division by zero FPU exceptions.
   */
  static void maskFpuDivZero();
  /**
   * Switch off signalling denormalized argument FPU exceptions.
   */
  static void maskFpuDenorm();
  /**
   * Switch off signalling invalid operation FPU exceptions.
   */
  static void maskFpuInvalid();
  //@}

};

}

#endif /* ThePEG_Debug_H */
