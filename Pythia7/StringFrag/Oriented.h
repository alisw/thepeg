// -*- C++ -*-
#ifndef PYTHIA7_Oriented_H
#define PYTHIA7_Oriented_H
// This is the declaration of the Oriented class.



//#include"FragConfig.h"

namespace Pythia7{

/**
 * The Oriented class mainly implements a static member theSide that
 * describes, in the fragmentation procedure, from which side of the
 * string a step is taken. The method to pick at random a side is
 * protected and only the LundFragHandler can change it. All other
 * classes can access to the current side in a fragmentation step by
 * the static <code>Dir()</code> method.
 *
 * In the future this class should be replaced by the standard
 * ThePEG::Direction class.
 */
class Oriented {

public:

  /** Enumeration of the directions */
  enum Side { left = -1, /**< Leftward direction. */
	      right = 1  /**< Rightward direction. */
  };

  /**
   * The current direction.
   */
  inline static int Dir();

  /**
   * The opposite of the current direction.
   */
  inline static int OppDir();

protected:

  /**
   * Pick the rightward (leftward) direction with probability \a rnd
   * (1-\a rnd).
   */
  inline static void pickSide(double rnd);

private:

  /** The current direction. */
  static Side  theSide; 

};

}
#include "Oriented.icc"

#endif // PYTHIA7_Oriented_H
