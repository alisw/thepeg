// -*- C++ -*-
#ifndef THEP8I_RndmEngine_H
#define THEP8I_RndmEngine_H
//
// This is the declaration of the RndmEngine class.
//

#include "TheP8I.h"
#include "Basics.h"

namespace TheP8I {

using namespace ThePEG;

/**
 * RndmEngine class inherits from the Pythia8 one and is used to allow
 * Pythia8 to use the same random number engine as ThePEG..
 */
class RndmEngine: public Pythia8::RndmEngine {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  RndmEngine();

  /**
   * The destructor.
   */
  virtual ~RndmEngine();
  //@}

public:

  /**
   * Generate a flat random number in the interval ]0,1[.
   */
  virtual double flat();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RndmEngine & operator=(const RndmEngine &);

};

}


#endif /* THEP8I_RndmEngine_H */
