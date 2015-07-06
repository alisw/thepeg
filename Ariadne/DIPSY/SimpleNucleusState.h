// -*- C++ -*-
#ifndef DIPSY_SimpleNucleusState_H
#define DIPSY_SimpleNucleusState_H
//
// This is the declaration of the SimpleNucleusState class.
//

#include "DipoleState.h"

namespace DIPSY {

using namespace ThePEG;

class SimpleNucleus;

/**
 * Here is the documentation of the SimpleNucleusState class.
 */
class SimpleNucleusState: public DipoleState {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline SimpleNucleusState() {}

  /**
   * The standard constructor.
   */
  SimpleNucleusState(const DipoleEventHandler & eh, SimpleNucleus & WF,
		     Energy plus, Energy minus, const vector<Point> & positions,
		     WFInfoPtr wfi);

  /**
   * The copy constructor.
   */
  inline SimpleNucleusState(const SimpleNucleusState & x) : DipoleState(x) {}

  /**
   * The destructor.
   */
  virtual ~SimpleNucleusState();
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SimpleNucleusState & operator=(const SimpleNucleusState &);

};

}

#endif /* DIPSY_SimpleNucleusState_H */
