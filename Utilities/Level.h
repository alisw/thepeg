// -*- C++ -*-
//
// Level.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Level_H
#define ThePEG_Level_H
// This is the declaration of the Level class.

#include "ThePEG/Config/ThePEG.h"

namespace ThePEG {

/**
 * Level is used to increment temporarily a given integer
 * variable. Everytime a Level object is created with a given integer
 * variable as argument, the variable will be incremented. When the
 * corresponding Level object is destroyed, the associated integer
 * variable is decremented again.
 *
 * @see HoldFlag
 */
template <typename T = int>
class Level {

public:

  /** Constructor taking an integer variable which is incremented. A
   *  reference to the variable will be stored. */
  Level(T & newLevel) : theLevel(++newLevel) {}

  /** Destructor decrementing the associated integer variable. */
  ~Level() { --theLevel; }

private:

  /** A reference to the integer variable to be decremmmented when
   *  this object is destroyed. */
  T & theLevel;

  /**
   * Default constructor is private and not implemented.
   */
  Level();

  /**
   * Copy constructor is private and not implemented.
   */
  Level(const Level &);

  /**
   * Assignment is private and not implemented.
   */
  Level & operator=(const Level &) = delete;

};

}

#endif /* ThePEG_Level_H */
